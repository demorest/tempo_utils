#! /usr/bin/env python
import re, struct, os
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
        elif len(line) > 80: 
            return "Tempo2"
        elif re.match("[0-9a-z@] ",line):
            return "Princeton"
        elif re.match("\S\S",line) and line[14]=='.':
            return "ITOA"
        elif re.match("  ",line) and line[41]=='.':
            return "Parkes"
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
        elif self.format=="Parkes" or self.format=="ITOA":
            raise RuntimeError, \
                "TOA format '%s' not implemented yet" % self.format

    def is_toa(self):
        """Return True if this is a valid TOA (ie not command/comment/etc)"""
        if self.format in ("Tempo2", "Princeton", "Parkes", "ITOA"):
            return True
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

class residual:
    def __init__(self,raw=None):
        self.mjd_bary = None
        self.res_phase = None
        self.res_us = None
        self.ophase = None
        self.rf_bary = None
        self.weight = None
        self.err_us = None
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
                self.rf_bary, self.weight, self.err_us, self.prefit_us,
                self.ddm, s2) = struct.unpack("=idddddddddi", raw)
        if s1 != 72 or s2 != 72:
            raise RuntimeError, "Error parsing residual block (s1=%d s2=%d)" \
                    % (s1, s2)
        # Change units:
        self.res_us *= 1e6
        self.prefit_us *= 1e6



def read_resid2_file(filename="resid2.tmp"):
    """Reads a tempo 'resid2.tmp' file and returns the result as a
    list of residual objects.  Each residual object has fields corresponding
    to the info in resid2.tmp:
      resid.mjd_bary    Barycentric MJD
      resid.res_phase   Residual in turns
      resid.res_us      Residual in us
      resid.ophase      Orbital phase in turns
      resid.rf_bary     Barycentric freq, MHz
      resid.weight      Weight of point in fit
      resid.err_us      TOA uncertainty, us
      resid.prefit_us   Prefit residual, us
      resid.ddm         DM correction from TOA line
    """
    f = open(filename, "r")
    resids = []
    r = f.read(80)
    while len(r)==80:
        resids += [residual(r)]
        r = f.read(80)
    f.close()
    return resids

class toalist(list):
    def get_resids(self,units='us'):
        if units=='us':
            return numpy.array([t.res.res_us for t in self if t.is_toa()])
        elif units=='phase':
            return numpy.array([t.res.res_phase for t in self if t.is_toa()])
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
        emax=0.0, ignore_EMAX=False):
    """Read the given filename and return a list of toa objects 
    parsed from it.  Will recurse to process INCLUDE-d files unless
    process_includes is set to False.  top is used internally for
    processing INCLUDEs and should always be set to True.  TOAs that are not
    marked by an in-file EMAX are filtered by emax > 0.0.  ignore_EMAX
    disables all in-file EMAXs."""
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
                if newtoa.command=="SKIP":
                    skip = True
                elif newtoa.command=="NOSKIP":
                    skip = False
                    continue
                if skip: continue
                if newtoa.command=="EMAX" and not ignore_EMAX:
                    emax = float(newtoa.arg)
                    continue
                if newtoa.command=="INFO":
                    read_toa_file.info = newtoa.arg
                if newtoa.command=="INCLUDE" and process_includes:
                    toas += read_toa_file(newtoa.arg,
                            ignore_blanks=ignore_blanks,top=False,emax=emax,
                            ignore_EMAX=ignore_EMAX)
                else:
                    toas += [newtoa]
            elif skip:
                continue
            elif newtoa.format in ("Blank", "Unknown"):
                if not ignore_blanks:
                    toas += [newtoa]
            else:
                newtoa.info = read_toa_file.info
                if emax>0.0 and newtoa.error>emax:
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
    if append:
        f = open(filename, "a")
    else:
        f = open(filename, "w")
    for t in toas:
        f.write(t.line + '\n')
    f.close()

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
def run_tempo(toas, parfile, show_output=False, get_output_par=False):
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
            print "tempo_utils.run_tempo: Adding 'MODE 1'"
            extra_cmds.insert(0,toa('MODE 1'))
        # Check if there are tempo2 TOAs but no FORMAT line
        if any([t.format=='Tempo2' for t in toas]) \
                and not any([t.command=='FORMAT' for t in toas]):
            print "tempo_utils.run_tempo: Adding 'FORMAT 1'"
            extra_cmds.insert(0,toa('FORMAT 1'))
        write_toa_file("pulsar.toa", extra_cmds + toas)
        cmd = "tempo -f pulsar.par pulsar.toa"
        if show_output==False:
            cmd += " > /dev/null"
        os.system(cmd)
        resids = read_resid2_file()
        toa_resid_match(toas, resids)
        lis = open("tempo.lis",'r').readlines()
        chi2 = float(lis[-1][14:23])
        ndof = float(lis[-1][24:30])
        rms = float((lis[-2].split())[4])
        # Capture output par file if needed
        if get_output_par:
            if psrname is not None:
                outparlines = open('%s.par' % psrname).readlines()
            else:
                outparlines = None
    finally:
        os.chdir(orig_dir)
    os.system("rm -rf %s" % temp_dir)
    if get_output_par:
        return (chi2,ndof,rms,outparlines)
    else:
        return (chi2,ndof,rms)

