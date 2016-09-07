#! /usr/bin/env python
import re, struct, os, string, shutil
import sys
import numpy

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
            raise RuntimeError, \
                "TOA format '%s' not implemented yet" % self.format

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
        raise RuntimeError,\
                "Record length does not match (s1=%d s2=%d len=%d)" % (ss[0],
                        ss[-1], dlen)
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
            raise RuntimeError, "Invalid raw residual block (len=%d)" \
                    % len(raw)
        (s1, self.mjd_bary, self.res_phase, self.res_us, self.ophase,
                self.rf_bary, self.weight, self.err_us, self.prefit_phase,
                self.ddm, s2) = struct.unpack("=idddddddddi", raw)
        if s1 != 72 or s2 != 72:
            raise RuntimeError, "Error parsing residual block (s1=%d s2=%d)" \
                    % (s1, s2)
        # Change units:
        self.res_us *= 1e6
        self.prefit_us = self.res_us / self.res_phase * self.prefit_phase


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
    f = open(filename, "r")
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
    f = open(filename, "r")
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

def toa_resid_match(toas, resids):
    """Match a list of residuals to a list of TOAs.  The res entry for each
    TOA is filled in appropriately.  resids list is altered along the way.
    Current algorithm is fairly dumb.."""
    toa_count = 0
    for toa in toas:
        if toa.is_toa():
            toa_count += 1
    if toa_count != len(resids):
        raise RuntimeError, "TOA count (%d) != residual count (%d)" \
                % (toa_count, len(resids))
    for toa in toas:
        if not toa.is_toa():
            continue
        toa.res = resids.pop(0)

import os, tempfile
def run_tempo(toas, parfile, show_output=False,
        get_output_par=False, gls=False, other_options=False,
        quiet=False,dcovfile=False):
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
                print "tempo_utils.run_tempo: Adding 'MODE 1'"
            extra_cmds.insert(0,toa('MODE 1'))
        # Check if there are tempo2 TOAs but no FORMAT line
        if any([t.format=='Tempo2' for t in toas]) \
                and not any([t.command=='FORMAT' for t in toas]):
            if not quiet:
                print "tempo_utils.run_tempo: Adding 'FORMAT 1'"
            extra_cmds.insert(0,toa('FORMAT 1'))
        write_toa_file("pulsar.toa", toalist(extra_cmds+toas))
        tempo_args = other_options if other_options else ""
        if gls: tempo_args += " -G"
        created_dcovfile = False
        if dcovfile:
            if os.path.exists(dcovfile): # See if dcovfile already exists
                shutil.copy(dcovfile,os.path.basename(dcovfile))
                if not any([l.startswith('DCOVFILE') for l in lines]):
                    lines.append('DCOVFILE ' + os.path.basename(dcovfile) + '\n')
            else: # dcovfile doesn't exist, so we will have tempo make it.
                tempo_args += " -C"
                created_dcovfile = True
                if any([l.startswith('DCOVFILE') for l in lines]):
                    lines = [l for l in lines if not l.startswith("DCOVFILE")]
        cmd = "tempo " + tempo_args + " -f pulsar.par pulsar.toa"
        if show_output==False:
            cmd += " > /dev/null"
        os.system(cmd)
        resids = read_resid2_file()
        toa_resid_match(toas, resids)
        lis = open("tempo.lis",'r').readlines()
        chi2_str = lis[-1][14:23]
        try:
            chi2 = float(chi2_str)
        except ValueError:
            chi2 = None
        ndof = float(lis[-1][24:30])
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
    finally:
        os.chdir(orig_dir)
    os.system("rm -rf %s" % temp_dir)
    if get_output_par:
        return (chi2,ndof,rms,outparlines)
    else:
        return (chi2,ndof,rms)

class polycos(list):
    """Collection of polyco sets."""

    def __init__(self,fname='polyco.dat'):
        fin = open(fname,'r')
        more = True
        while more:
            try:
                self.append(polyco(fin=fin))
            except:
                more = False

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

    def phase_and_freq(self,mjd,fmjd=0.0):
        dt_min = (float(mjd - self.imjd) + (fmjd - self.fmjd))*1440.0
        if abs(dt_min) > self.span/2.0:
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
