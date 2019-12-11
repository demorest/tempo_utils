#! /usr/bin/env python
import re, struct, os, string, shutil
import subprocess
import sys
import numpy
import tempfile
from collections import namedtuple

toa_commands = ("DITHER", "EFAC", "EMAX", "EMAP", "EMIN", "EQUAD", "FMAX",
        "FMIN", "INCLUDE", "INFO", "JUMP", "MODE", "NOSKIP", "PHA1", "PHA2",
        "PHASE", "SEARCH", "SIGMA", "SIM", "SKIP", "TIME", "TRACK", "ZAWGT",
        "FORMAT", "EFLOOR")

def toa_format(line):
    """Identifies a TOA line as one of the following types:  Comment, Command,
    Blank, Tempo2, Princeton, ITOA, Parkes, Unknown."""
    try:
        if line[0]=='C' or line[0]=='#':
            return "Comment"
        elif line.startswith(toa_commands):
            return "Command"
        elif re.match("^\s+$", line):
            return "Blank"
        elif re.match("[0-9a-z@] ",line):
            return "Princeton"
        elif (re.match("\S\S",line) and line[14]=='.'
                and re.match("[a-zA-Z]{2}",line[57:59])):
            return "ITOA"
        elif re.match(" ",line) and line[41]=='.':
            return "Parkes"
        elif len(line.split())>=5:
            return "Tempo2"
        else:
            return "Unknown"
    except:
        return "Unknown"

class toa:
    def __init__(self, line=None):
        self.line = line.rstrip('\n')
        if self.line != None:
            self.parse_line()

    def __repr__(self):
        return self.line

    def parse_line(self):
        self.command = None
        self.arg = None
        self.site = None
        self.freq = None
        self.mjd = None
        self.imjd = None
        self.fmjd = None
        self.error = None
        self.ddm = None
        self.res = None
        self.fname = None
        self.info = None
        self.flags = {}
        self.format = toa_format(self.line)
        if self.format=="Command":
            tmp = self.line.split()
            if len(tmp)==1:
                self.command = tmp[0]
            elif len(tmp)>1:
                self.command = tmp[0]
                self.arg = tmp[1]
        elif self.format=="Princeton":
            self.site = self.line[0]
            self.freq = float(self.line[15:24])
            self.mjd = float(self.line[24:44])
            self.imjd = int(self.mjd)
            if self.line[29]=='.':
                self.fmjd = float(self.line[29:44])
            elif self.line[30]=='.':
                self.fmjd = float(self.line[30:44])
            self.error = float(self.line[44:53])
            try:
                self.ddm = float(self.line[68:78])
            except:
                self.ddm = 0.0
        elif self.format=="Parkes":
            self.freq = float(self.line[25:34])
            self.mjd = float(self.line[34:55])
            self.imjd = int(self.mjd)
            self.fmjd = float(self.line[41:55])
            # TODO phase offset
            self.error = float(self.line[63:71])
            self.site = self.line[79]
        elif self.format=="Tempo2":
            # This could use more error catching...
            fields = self.line.split()
            self.fname = fields.pop(0)
            self.freq = float(fields.pop(0))
            mjdstr = fields.pop(0)
            self.mjd = float(mjdstr)
            self.imjd = int(self.mjd)
            self.fmjd = float(mjdstr[mjdstr.find('.'):])
            self.error = float(fields.pop(0))
            self.site = fields.pop(0)
            # All the rest should be flags
            for i in range(0,len(fields),2):
                self.flags[fields[i].lstrip('-')] = fields[i+1]
        elif self.format=="ITOA":
            raise RuntimeError( 
                    "TOA format '%s' not implemented yet" % self.format)

    def is_toa(self):
        """Return True if this is a valid TOA (ie not command/comment/etc)"""
        if self.format in ("Tempo2", "Princeton", "Parkes", "ITOA"):
            return True
        else:
            return False

    def is_commented_toa(self):
        """Return True if this is a commented out TOA line"""
        if self.format == "Comment":
            try:
                return toa(self.line.lstrip('C# ')).is_toa()
            except:
                return False
        else:
            return False

    def flag(self,flagname):
        """Return the value of the flag, or None if flag not present."""
        try:
            rv = self.flags[flagname]
        except KeyError:
            rv = None
        return rv

    def comment(self):
        """Comment out this TOA line."""
        if self.format != 'Comment':
            self.__init__('C ' + self.line)

    def uncomment(self):
        """Uncomment this TOA line."""
        if self.format == 'Comment':
            self.__init__(self.line.lstrip('C# '))

def _unpack_record(raw,fmt=None):
    dlen = len(raw) - 8
    if fmt is None:
        fmt = 'd' * (dlen/8)
    ss = struct.unpack('=i'+fmt+'i',raw)
    if ss[0]!=dlen or ss[-1]!=dlen:
        raise RuntimeError( 
                "Record length does not match (s1=%d s2=%d len=%d)" % (ss[0],
                        ss[-1], dlen))
    return ss[1:-1]

class residual:
    def __init__(self,raw=None):
        self.mjd_bary = None
        self.res_phase = None
        self.res_us = None
        self.ophase = None
        self.rf_bary = None
        self.weight = None
        self.err_us = None
        self.prefit_phase = None
        self.prefit_us = None
        self.ddm = None
        if raw!=None:
            self.parse_raw(raw)

    def __repr__(self):
        if self.mjd_bary==None:
            return "Uninitialized residual"
        return "Residual at MJD %.10f, %.3f MHz: %+.3f +/- %.3f us" \
                % (self.mjd_bary, self.rf_bary, self.res_us, self.err_us)

    def parse_raw(self,raw):
        if len(raw) != 80:
            raise RuntimeError("Invalid raw residual block (len=%d)" 
                    % len(raw))
        (s1, self.mjd_bary, self.res_phase, self.res_us, self.ophase,
                self.rf_bary, self.weight, self.err_us, self.prefit_phase,
                self.ddm, s2) = struct.unpack("=idddddddddi", raw)
        if s1 != 72 or s2 != 72:
            raise RuntimeError("Error parsing residual block (s1=%d s2=%d)" 
                    % (s1, s2))
        # Change units:
        self.res_us *= 1e6
        self.prefit_us = (self.res_us / self.res_phase * self.prefit_phase 
                if self.res_phase!=0.0 else 0.0)


def read_resid2_file(filename="resid2.tmp"):
    """Reads a tempo 'resid2.tmp' file and returns the result as a
    list of residual objects.  Each residual object has fields corresponding
    to the info in resid2.tmp:
      resid.mjd_bary      Barycentric MJD
      resid.res_phase     Residual in turns
      resid.res_us        Residual in us
      resid.ophase        Orbital phase in turns
      resid.rf_bary       Barycentric freq, MHz
      resid.weight        Weight of point in fit
      resid.err_us        TOA uncertainty, us
      resid.prefit_phase  Prefit residual in turns
      resid.prefit_us     Prefit residual, us
      resid.ddm           DM correction from TOA line
    """
    f = open(filename, "rb")
    resids = []
    r = f.read(80)
    while len(r)==80:
        resids += [residual(r)]
        r = f.read(80)
    f.close()
    return resids

def read_design_matrix(filename='design.tmp'):
    """Reads a tempo 'design.tmp' file and returns the result as a
    numpy array, dims (ntoa,nparam)"""
    f = open(filename, "rb")
    # First record should be two ints giving array size
    (ntoa,nparam) = _unpack_record(f.read(16),'ii')
    result = numpy.zeros((ntoa,nparam+1))
    # First two elements are time and weight, which we will ignore
    # here.
    rsize = (nparam+2)*8 + 8
    for i in range(ntoa):
        result[i,:-1] = _unpack_record(f.read(rsize))[2:]
    # Constant phase term
    result[:,-1] = 1.0
    f.close()
    return result

class toalist(list):

    def write(self, filename, append=False):
        if append:
            f = open(filename, "a")
        else:
            f = open(filename, "w")
        for t in self:
            f.write(t.line + '\n')
        f.close()

    def get_ntoa(self,commented=False):
        ntoa = sum(t.is_toa() for t in self)
        if commented: ntoa += sum(t.is_commented_toa() for t in self)
        return ntoa
    def get_resids(self,units='us'):
        if units=='us':
            return numpy.array([t.res.res_us for t in self if t.is_toa()])
        elif units=='phase':
            return numpy.array([t.res.res_phase for t in self if t.is_toa()])
    def get_prefit(self,units='us'):
        if units=='us':
            return numpy.array([t.res.prefit_us for t in self if t.is_toa()])
        elif units=='phase':
            return numpy.array([t.res.prefit_phase for t in self if t.is_toa()])
    def get_orb_phase(self):
        return numpy.array([t.res.ophase for t in self if t.is_toa()])
    def get_resid_err(self):
        return numpy.array([t.res.err_us for t in self if t.is_toa()])
    def get_freq(self):
        return numpy.array([t.freq for t in self if t.is_toa()])
    def get_mjd(self):
        return numpy.array([t.mjd for t in self if t.is_toa()])
    def get_err(self):
        return numpy.array([t.error for t in self if t.is_toa()])
    def get_flag(self,flag,f=lambda x: x):
        return numpy.array([f(t.flags[flag]) for t in self if t.is_toa()])
    def get_chi2(self):
        # NOTE: this will not take EQUAD/EFAC/etc into account
        #x = self.get_resids()/self.get_err()
        # This should:
        x = self.get_resids()/self.get_resid_err()
        return (x**2).sum()
    def get_group_chi2(self,flag,flagval,reduced=False,remove_mean=False):
        r = numpy.array([t.res.res_us for t in self if t.is_toa() and t.flags[flag]==flagval])
        e = numpy.array([t.res.err_us for t in self if t.is_toa() and t.flags[flag]==flagval])
        if remove_mean:
            w = 1.0/(e*e)
            rm = (r*w).sum() / w.sum()
            r = r - rm
        x = r / e
        if reduced:
            return (x**2).mean()
        else:
            return (x**2).sum()

def read_toa_file(filename,process_includes=True,ignore_blanks=True,top=True,
        emax=0.0, process_emax=True, fmin=0.0, process_fmin=True,
        fmax=0.0, process_fmax=True, process_skips=True, convert_skips=False):
    """Read the given filename and return a list of toa objects
    parsed from it.   Options:

        process_includes: if True, read TOA lines from any INCLUDED files.

        ignore_blanks: if True, blank TOA lines are ignored.

        top: Used internally for recursion, do not use this.

        emax: Set to apply a certain EMAX value as TOAs are read.

        fmin: Set to apply a certain FMIN value as TOAs are read.

        fmax: Set to apply a certain FMAX value as TOAs are read.

        process_emax: if True, EMAX statements are applied as TOAs are read.

        process_fmin: if True, FMIN statements are applied as TOAs are read.

        process_fmax: if True, FMAX statements are applied as TOAs are read.

        process_skips: if True, SKIP'd sections are not read in.

        convert_skips: if True, SKIP'd lines are converted to comments.
    """
    # Change to directory where top-level file is in order to process
    # relative include paths correctly.  This might break if there
    # are multiple levels of include...
    orig_dir = None
    try:
        if top and ('/' in filename):
            orig_dir = os.getcwd()
            work_dir = '/'.join(filename.split('/')[:-1])
            filename = filename.split('/')[-1]
            os.chdir(work_dir)
        f = open(filename, "r")
        toas = toalist([])
        skip = False
        if top: read_toa_file.info = None
        for l in f.readlines():
            newtoa = toa(l)
            if newtoa.format=="Command":
                if newtoa.command=="SKIP" and process_skips:
                    skip = True
                elif newtoa.command=="NOSKIP" and process_skips:
                    skip = False
                    continue
                if skip:
                    if convert_skips: newtoa.comment()
                    else: continue
                if newtoa.command=="EMAX" and process_emax:
                    emax = float(newtoa.arg)
                    continue
                if newtoa.command=="FMIN" and process_fmin:
                    fmin = float(newtoa.arg)
                    continue
                if newtoa.command=="FMAX" and process_fmax:
                    fmax = float(newtoa.arg)
                    continue
                if newtoa.command=="INFO":
                    read_toa_file.info = newtoa.arg
                if newtoa.command=="INCLUDE" and process_includes:
                    toas += read_toa_file(newtoa.arg,
                            ignore_blanks=ignore_blanks,top=False,
                            emax=emax, process_emax=process_emax,
                            fmin=fmin, process_fmin=process_fmin,
                            fmax=fmax, process_fmax=process_fmax,
                            process_skips=process_skips)
                else:
                    toas += [newtoa]
            elif skip:
                if convert_skips: newtoa.comment()
                else: continue
            elif newtoa.format in ("Blank", "Unknown"):
                if not ignore_blanks:
                    toas += [newtoa]
            else:
                newtoa.info = read_toa_file.info
                if emax>0.0 and newtoa.error>emax:
                    pass
                elif fmin>0.0 and newtoa.freq<fmin:
                    pass
                elif fmax>0.0 and newtoa.freq>fmax:
                    pass
                else:
                    toas += [newtoa]
        f.close()
        if top and (orig_dir is not None):
            os.chdir(orig_dir)
        return toas
    except:
        # If something went wrong, make sure we go back to the original
        # directory, then re-raise the exception
        if orig_dir is not None:
            os.chdir(orig_dir)
        raise

def write_toa_file(filename, toas, append=False):
    toas.write(filename,append=append)

def toa_match(toa1, toa2, freq_tol=0.1, time_tol=0.001):
    """Return true if two toa objects match within the given tolerances.
    Always returns false for command/etc lines.  freq_tol is in MHz and
    time_tol is in days."""
    if not toa1.is_toa(): return False
    if not toa2.is_toa(): return False
    if abs(toa2.freq - toa1.freq) > freq_tol: return False
    # Don't need super-high accuracy for the date comp
    if abs(toa2.mjd - toa1.mjd) > time_tol: return False
    return True

def toa_list_match(toas1, toas2, freq_tol=0.1, time_tol=0.001):
    """Match up two sets of TOAs.  Returns two lists of the same length
    containing only the TOAs that match."""
    out1 = toalist([])
    out2 = toalist([])
    for toa1 in toas1:
        for toa2 in toas2:
            if toa_match(toa1, toa2, freq_tol, time_tol):
                out1 += [toa1]
                out2 += [toa2]
                #toas1.remove(toa1)
                #toas2.remove(toa2)
                break
    return (out1, out2)

def toa_resid_match(toas, resids, phi=None):
    """Match a list of residuals to a list of TOAs.  The res entry for each
    TOA is filled in appropriately.  resids list is altered along the way.
    Current algorithm is fairly dumb.."""
    toa_count = 0
    for toa in toas:
        if toa.is_toa():
            toa_count += 1
    if toa_count != len(resids):
        raise RuntimeError("TOA count (%d) != residual count (%d)" 
                % (toa_count, len(resids)))
    for toa in toas:
        if not toa.is_toa():
            continue
        toa.res = resids.pop(0)
        if phi is not None:
            toa.res.phi = phi.pop(0)

import os, tempfile
def run_tempo(toas, parfile, show_output=False,
        get_output_par=False, gls=False, other_options=False,
        quiet=False, dcovfile=False, get_phisun=False, inabspulsenum=False,
        matrixfile=False):
    """Run tempo on the given TOA list using the given parfile.  Residuals
    are read and filled into the toa structs on successful completion."""
    orig_dir = os.getcwd()
    try:
        temp_dir = tempfile.mkdtemp(prefix="tempo")
        psrname = None
        try:
            lines = open(parfile,'r').readlines()
        except:
            lines = parfile
        for l in lines:
            if l.startswith('PSR'):
                psrname = l.split()[1]
        #os.system("cp %s %s/pulsar.par" % (parfile, temp_dir))
        open("%s/pulsar.par" % temp_dir,'w').writelines(lines)
        os.chdir(temp_dir)
        extra_cmds = toalist([])
        # Always add mode 1 if it's not there
        if not any([t.command=='MODE' for t in toas]):
            if not quiet:
                print("tempo_utils.run_tempo: Adding 'MODE 1'")
            extra_cmds.insert(0,toa('MODE 1'))
        # Check if there are tempo2 TOAs but no FORMAT line
        if any([t.format=='Tempo2' for t in toas]) \
                and not any([t.command=='FORMAT' for t in toas]):
            if not quiet:
                print("tempo_utils.run_tempo: Adding 'FORMAT 1'")
            extra_cmds.insert(0,toa('FORMAT 1'))
        write_toa_file("pulsar.toa", toalist(extra_cmds+toas))
        tempo_args = other_options if other_options else ""
        if gls: tempo_args += " -G"
        if get_phisun: tempo_args += " -a"
        created_dcovfile = False
        if dcovfile:
            if os.path.exists(dcovfile): # See if dcovfile already exists
                shutil.copy(dcovfile,os.path.basename(dcovfile))
                if not any([l.startswith('DCOVFILE') for l in lines]):
                    open("pulsar.par",'a').writelines('DCOVFILE ' + os.path.basename(dcovfile) + '\n')
            else: # dcovfile doesn't exist, so have tempo make it.
                tempo_args += " -C"
                created_dcovfile = True
        if inabspulsenum:
            if os.path.exists(inabspulsenum):
                shutil.copy(inabspulsenum,os.path.basename(inabspulsenum))
                tempo_args += " -ni " + os.path.basename(inabspulsenum)
            else:
                print(inabspulsenum," does not exist")
        cmd = "tempo " + tempo_args + " -f pulsar.par pulsar.toa"
        if show_output==False:
            cmd += " > /dev/null"
        os.system(cmd)
        resids = read_resid2_file()
        if get_phisun: 
            phisun = list(numpy.loadtxt('phisun.tmp'))
        else:
            phisun = None
        toa_resid_match(toas, resids, phisun)
        lis = open("tempo.lis",'r').readlines()
        idx1 = lis[-1].index(':')+1
        idx2 = lis[-1].index('=')
        chi2_str = lis[-1][idx1:idx2].split('/')
        try:
            chi2 = float(chi2_str[0])
        except ValueError:
            chi2 = None
        ndof = float(chi2_str[1])
        rms_str = lis[-2][63:74]
        try:
            rms = float(rms_str)
        except ValueError:
            rms = None
        # Capture output par file if needed
        if get_output_par:
            if psrname is not None:
                outparlines = open('%s.par' % psrname).readlines()
            else:
                outparlines = None
        if created_dcovfile:
            os.rename("datacov.tmp",dcovfile)
        if matrixfile:
            shutil.copy("matrix.tmp",os.path.join(orig_dir,matrixfile))
    finally:
        os.chdir(orig_dir)
    os.system("rm -rf %s" % temp_dir)
    if get_output_par:
        return (chi2,ndof,rms,outparlines)
    else:
        return (chi2,ndof,rms)

class polycos(list):
    """Collection of polyco sets."""

    def __init__(self,fname='polyco.dat',tempo_args=None,tempo_output=None):
        fin = open(fname,'r')
        more = True
        self.tempo_args = tempo_args
        self.tempo_output = tempo_output
        while more:
            try:
                self.append(polyco(fin=fin))
            except:
                more = False

    @classmethod
    def generate_from_polyco(cls, parfile, polyco0):
        """Generate a polyco set with the same parameters as an existing
        single polyco instance (polyco0)."""
        return cls.generate(parfile, 
                polyco0.site,
                polyco0.imjd + polyco0.fmjd - polyco0.span/2.0/1440.0,
                polyco0.span/1440.0, 
                polyco0.ncoeff,
                polyco0.obsfreq)

    @classmethod
    def generate(cls, parfile, site, mjd_start, tobs, 
            ncoeff=15, freq=1420.0, outfile=None):
        """Generate a polyco set from parfile.  If outfile is not none,
        the polycos will be saved in the given file name, otherwise
        deleted after they are read."""
        # See if parfile is lines or a path
        try:
            parlines = open(parfile,'r').readlines()
        except:
            parlines = parfile
        for l in parlines:
            if l.startswith('PSR'):
                psrname = (l.split()[1]).lstrip('BJ')
        tmpdir = None
        origdir = os.getcwd()
        try:
            tmpdir = tempfile.mkdtemp(prefix='polyco')
            open("%s/pulsar.par" % tmpdir,'w').writelines(parlines)
            if outfile is not None: fullout=os.path.abspath(outfile)
            os.chdir(tmpdir)
            try:
                tempo_args = ['tempo', '-f', 'pulsar.par',
                        '-ZPSR=%s' % psrname,
                        '-ZNCOEFF=%d' % ncoeff,
                        '-ZFREQ=%.6f' % freq,
                        '-ZSITE=%s' % site,
                        '-ZSTART=%.8f' % mjd_start,
                        '-ZTOBS=%.4f' % tobs]
                tempo_output = subprocess.check_output(tempo_args,
                        stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                tempo_output = e.output
            if outfile is not None:
                shutil.copy('polyco.dat', fullout)
            return cls(tempo_args=tempo_args,tempo_output=tempo_output)
        finally:
            os.chdir(origdir)
            shutil.rmtree(tmpdir)

    def match(self,mjd,fmjd=0.0):
        for p in self:
            try:
                foo = p.phase_and_freq(mjd,fmjd)
                return p
            except RuntimeError:
                pass
        raise RuntimeError('No matching polycos found (mjd=%.10f)'%(mjd+fmjd))

    def phase_and_freq(self,mjd,fmjd=0.0):
        for p in self:
            try:
                return p.phase_and_freq(mjd,fmjd)
            except RuntimeError:
                pass
        raise RuntimeError('No matching polycos found (mjd=%.10f)'%(mjd+fmjd))

    def phase(self,mjd,fmjd=0.0):
        (phase,freq) = self.phase_and_freq(mjd,fmjd)
        return phase

    def freq(self,mjd,fmjd=0.0):
        (phase,freq) = self.phase_and_freq(mjd,fmjd)
        return freq

class polyco:
    """Single polyco set."""

    def __init__(self, fname='polyco.dat', fin=None):
        # TODO: error checking
        if fin is None:
            if fname is None:
                # Init everything with dummy values
                self.src = 'psr'
                self.imjd = 0
                self.fmjd = 0.0
                self.dm = 0.0
                self.earth_z4 = 0.0
                self.log10_rms = 10.0
                self.rphase = 0.0
                self.rfreq = 0.0
                self.site = '0'
                self.span = 0.0
                self.ncoeff = 1
                self.obsfreq = 0.0
                self.ophase = None
                self.coeffs = [0.0,]
                return
            else:
                fin = open(fname,'r')
        line1 = fin.readline().strip().split()
        line2 = fin.readline().strip().split()
        self.src = line1[0]
        (s_imjd, s_fmjd) = line1[3].split('.')
        self.imjd = int(s_imjd)
        self.fmjd = float('.' + s_fmjd)
        self.dm = float(line1[4])
        self.earth_z4 = float(line1[5])
        self.log10_rms = float(line1[6])
        self.rphase = float(line2[0])
        self.rfreq = float(line2[1])
        self.site = line2[2]
        self.span = float(line2[3])
        self.ncoeff = int(line2[4])
        self.obsfreq = float(line2[5])
        try:
            self.ophase = float(line2[6])
        except IndexError:
            self.ophase = None
        nclines = self.ncoeff / 3
        if self.ncoeff % 3:
            nclines += 1
        coeff_str = ''
        for i in range(nclines):
            coeff_str += fin.readline().strip() + ' '
        self.coeffs = map(float,coeff_str.replace('D','e').split())

    def as_string(self):
        # Return the info formatted as a polyco block
        hh = int(self.fmjd*24.0)
        mm = int(self.fmjd*24.0*60.0 - hh*60.0)
        ss = self.fmjd*86400.0 - hh*3600.0 - mm*60.0
        mjd_str = str(self.imjd) + ('%.11f'%self.fmjd).lstrip('0')
        line1 = '%-10s %9s%11.2f%20s%21.6f %6.3f%7.3f' % (
                self.src,
                'DD-MMM-YY',
                hh*10000 + mm*100 + ss,
                mjd_str,
                self.dm,
                self.earth_z4,
                self.log10_rms
                )
        line2 = '%20.6f%18.12f%5s%5d%5d%10.3f%7.4f' % (
                self.rphase,
                self.rfreq,
                self.site,
                int(self.span),
                self.ncoeff,
                self.obsfreq,
                self.ophase if self.ophase is not None else 0.0
                )
        lines = [line1, line2]
        for i in range(0,self.ncoeff,3):
            cline = ''
            for c in self.coeffs[i:i+3]:
                cline += ('%25.17E' % c).replace('E','D')
            lines += [cline,]
        return string.join(lines,'\n') + '\n'

    def phase_and_freq(self,mjd,fmjd=0.0):
        dt_min = (float(mjd - self.imjd) + (fmjd - self.fmjd))*1440.0
        if abs(dt_min) > (1.01*self.span/2.0):
            raise RuntimeError('MJD outside polyco span (dt=%.1f min)' % dt_min)
        freq = 0.0
        phase = self.coeffs[self.ncoeff-1]
        for i in range(self.ncoeff-1,0,-1):
            phase = dt_min*phase + self.coeffs[i-1]
            freq = dt_min*freq + float(i)*self.coeffs[i]
        freq = self.rfreq + freq/60.0
        phase += self.rphase + dt_min*60.0*self.rfreq
        return (phase, freq)

    def phase(self,mjd,fmjd=0.0):
        (phase,freq) = self.phase_and_freq(mjd,fmjd)
        return phase

    def freq(self,mjd,fmjd=0.0):
        (phase,freq) = self.phase_and_freq(mjd,fmjd)
        return freq

class parfile(object):

    # simple struct for DMX
    _dmx_range = namedtuple('dmx_range',
            ['idx','val','ep','r1','r2','f1','f2'])

    # Regex for matching T2-style JUMPs
    _jump_re = re.compile('(JUMP\s+\S+\s+\S+)\s')

    def __init__(self, fname):
        # Try it as a filename first
        try:
            self.lines = open(fname,'r').readlines()
            return
        except TypeError:
            pass
        # Try it as a file object next
        try:
            self.lines = fname.readlines()
            return
        except AttributeError:
            pass
        # Try it as another parfile object next (makes a copy)
        try:
            self.lines = list(fname.lines)
            return
        except:
            pass
        # Try it as a list of strings
        # This one might need more error checking...
        try:
            self.lines = list(fname)
        except:
            raise TypeError("parfile can't be constructed from type '%s'" % 
                    str(type(fname)))

    def write(self, fname):
        """Write out the par file to file fname."""
        open(fname,'w').writelines(self.lines)

    def run_fit(self, toas, gls=None):
        """Run tempo using the parfile and TOAs, return results as a new 
        parfile instance."""
        if gls is None:
            if 'ECORR' in self.keys or 'RNAMP' in self.keys:
                gls = True
            else:
                gls = False
        chi2, dof, rms, lines = run_tempo(toas, self.lines, 
                get_output_par=True, gls=gls)
        result = parfile(lines)
        result.chi2 = chi2
        result.dof = dof
        result.rms = rms
        return result

    @staticmethod
    def _is_comment(parline):
        if parline.startswith('#') or parline.startswith('C '):
            return True
        else:
            return False

    @staticmethod
    def _is_blank(parline):
        return re.match("^\s+$", parline)

    @staticmethod
    def _fortran_float(strval):
        return float(strval.replace('D','e'))

    @property
    def keys(self):
        keys = []
        for l in self.lines:
            if not self._is_comment(l) and not self._is_blank(l):
                jump = self._jump_re.match(l)
                if jump:
                    keys.append(jump.group(1))
                else:
                    keys.append(l.split()[0])
        return keys

    def _iline(self,key):
        if not key in self.keys:
            raise KeyError("key '%s' not found" % key)
        for i, l in enumerate(self.lines):
            if l.startswith(key+' '):
                return i

    def line(self, key, strip=False):
        ltmp = self.lines[self._iline(key)]
        if strip:
            ltmp = ltmp.replace(key,'',1)
            return ltmp.strip()
        return self.lines[self._iline(key)].rstrip()

    def val(self, key, dtype=str):
        l = self.line(key,strip=True)
        if dtype=='float' or dtype==float: 
            dtype_func = self._fortran_float
        else: 
            dtype_func = dtype
        return dtype_func(l.split()[0])

    def err(self, key):
        ltmp = self.line(key,strip=True)
        vals = ltmp.split()
        nval = len(vals)
        if nval>=3:
            return self._fortran_float(vals[2])
        else:
            return None

    def set_val(self,key,val):
        """Set a new value for this parameter."""
        # Not much checking of whether this is allowed for the given
        # param.  val should be a string
        idx = self._iline(key)
        ltmp = self.line(key,strip=True)
        vals = ltmp.split()
        vals[0] = val
        ltmp = string.join([key,] + vals,' ') + '\n'
        self.lines[idx] = ltmp

    def remove(self,key):
        """Remove this parameter from the parfile (delete line).  No error
        is raised if the parameter is not present."""
        try:
            self.lines.pop(self._iline(key))
        except KeyError:
            pass

    def add(self,line):
        """Add the line to the parfile."""
        self.lines += [line.rstrip() + '\n',]

    def no_fit(self):
        """Turn all fit params off.  Note, not fully general yet."""
        self.lines = [l.replace(' 1 ',' 0 ') for l in self.lines]

    def is_fit(self,key):
        """Return whether a given param is being fit."""
        ltmp = self.line(key,strip=True)
        vals = ltmp.split()
        nval = len(vals)
        if nval>=2:
            return vals[1]=='1'
        return False

    def set_fit(self,key,fit=True):
        """Turn the fit for this parameter on or off."""
        # Not much checking of whether this is allowed for the given
        # param.
        if fit and self.is_fit(key): return
        if not fit and not self.is_fit(key): return
        idx = self._iline(key)
        ltmp = self.line(key,strip=True)
        vals = ltmp.split()
        nval = len(vals)
        if nval==1 and fit:
            ltmp = string.join([key,]+vals+['1 \n',],' ')
        elif nval>=2:
            vals[1] = '1' if fit else '0'
            ltmp = string.join([key,]+vals,' ') + '\n'
        self.lines[idx] = ltmp

    def dmx(self,dmxidx):
        """Return all dmx info for the given DMX range."""
        val = self.val('DMX_%04d' % dmxidx, float)
        stuff = [dmxidx,val,]
        for k in 'DMXEP', 'DMXR1', 'DMXR2', 'DMXF1', 'DMXF2':
            try:
                stuff += [self.val('%s_%04d' % (k, dmxidx), float),]
            except KeyError:
                stuff += [None,]
        return self._dmx_range(*stuff)

    @property
    def dmx_indices(self):
        """Return a list of all dmx indices in the file."""
        return sorted([int(k.split('_')[1])
            for k in self.keys 
            if k.startswith('DMX_')])

