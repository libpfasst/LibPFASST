#!/usr/bin/python
# Copyright (C) 2006-2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import re, sys, os, getopt

def usage():
    print """f90_mod_deps.py [options] <file ...>

Generates dependencies for the given Fortran 90 source files based on
module and use statements in them. Options are:

  -h, --help     This help output.

  -o, --output <file>
                 Specify the output file. If unspecified then output
                 is to stdout.

  -d, --dep-re <regexp>
                 Regular expression to match against each module name
                 within use statements. Defaults to matching
                 everything.

  -D, --dep-template <template>
                 Template expression for the dependency produced from
                 each module matched by the dep-re regular expression.

  -m, --mod-re <regexp>
                 Regular expression to match against each module name
                 within module definition statements. Defaults to
                 matching everything.

  -M, --mod-template <template>
                 Template expression for the dependency target
                 produced from each module matched by the mod-re
                 regular expression.

  -op, --o-prefix <prefix>

  -v, --verbose  Turn on verbose debugging output.

For a discussion of managing Fortran 90 dependencies see:
http://tableau.stanford.edu/~mwest/group/Fortran_90_Module_Dependencies

Example:
f90_mod_deps.py --output src/myfile.deps --dep-re "(pmc_.*)" \\
      --dep-template "src/\1.mod" --mod-re "(.*)" \\
      --mod-template "src/\1.mod" src/myfile.f90
"""

# default options
class Opts:
    output = None
    dep_re = "(.*)"
    dep_template = "\\1.mod"
    mod_re = "(.*)"
    mod_template = "\\1.mod"
    o_prefix = ""
    verbose = False

def get_deps_and_mods(filename, opts):
    ''' some changes by ab (2010.11.30):
    below is the code used for testing:
    	
    USE treevars, ONLY: x, & ! comment for testing
       yy => y, &
       z
    use utils ! more comment for testing
    
    The FORTRAN 2008 standard (ftp://ftp.nag.co.uk/sc22wg5/N1801-N1850/N1830.pdf) allows only one module per use statement.
    After a "use modname" only a ", only:..." or something like ", x => y ..." can be used.
    The standard allows several ";"-seperated use-statements in one line, but since it is very hard to account for this and nobody uses it i have not implemented it.
    The standard also allows a "module-nature" before the module name, but again i have not implemented it
    
    '''

    if opts.verbose:
	sys.stderr.write("Processing %s\n" % filename)
    deps = []
    mods = []
    f = open(filename)
    if not f:
	print "ERROR: unable to open %s%s" % filename
	sys.exit(1)
    use_line_re = re.compile("^\s*use\s+([a-z][a-z0-9_]+)", re.IGNORECASE)				# 1: new by ab: takes only one valid module name ([a-z][a-z0-9_]+), almost following
    													#  the fortran 2008 standard (ftp://ftp.nag.co.uk/sc22wg5/N1801-N1850/N1830.pdf, p.273)
    													#  but ignoring the "module-nature" and the possibility to write several ";"-seperated
    													#  USE-statements in one line. now ignoring case.
    													#  this regexp is not checking for a valid use-statement, that's the compilers job.
    mod_line_re = re.compile("^\s*module\s+([a-z][a-z0-9_]+)\s*$", re.IGNORECASE)			# 3: new by ab: now matching only valid module names (see 1:). ignoring case
    split_re = re.compile("\s*,\s*")
    dep_re = re.compile(opts.dep_re)
    mod_re = re.compile(opts.mod_re)
    for line in f:
	line = re.sub("!.*","",line)									# 4: new by ab: first remove every comment
	match = use_line_re.search(line)
	if match:
	    item = match.group(1)
	    if item:
		if opts.verbose:
		    sys.stderr.write("use: %s\n" % item)
		match = dep_re.match(item)
		if match:
		    dep = match.expand(opts.dep_template)
		    if opts.verbose:
		        sys.stderr.write("matched to: %s\n" % dep)
		    if dep not in deps:
		        deps.append(dep)
	
	else:
	    match = mod_line_re.search(line)
	    if match:
		mod_name = match.group(1)
		if opts.verbose:
		    sys.stderr.write("module: %s\n" % mod_name)
		match = mod_re.match(mod_name)
		if match:
		    mod = match.expand(opts.mod_template)
		    if opts.verbose:
			sys.stderr.write("matched to: %s\n" % mod)
		    if mod not in mods:
			mods.append(mod)
    f.close()
    return (deps, mods)

def write_deps(outf, filename, deps, opts):
    filebase, fileext = os.path.splitext(filename)
    outf.write("%s%s.o: %s\n" % (opts.o_prefix, os.path.basename(filebase), " ".join(deps)))

def process_args():
    try:
	opts, args = getopt.getopt(sys.argv[1:], "ho:d:D:m:M:v",
				   ["help", "output=", "dep-re=",
				    "dep-template=", "mod-re=",
				    "mod-template=", "o-prefix=", "verbose"])
    except getopt.GetoptError:
	print "ERROR: invalid commandline options"
	usage()
	sys.exit(1)
    myopts = Opts()
    for o, a in opts:
	if o in ("-h", "--help"):
	    usage()
	    sys.exit()
	if o in ("-o", "--output"):
	    myopts.output = a
	if o in ("-d", "--dep-re"):
	    myopts.dep_re = a
	if o in ("-D", "--dep-template"):
	    myopts.dep_template = a
	if o in ("-m", "--mod-re"):
	    myopts.mod_re = a
	if o in ("-M", "--mod-template"):
	    myopts.mod_template = a
	if o in ("-op", "--o-prefix"):
	    myopts.o_prefix = a
	if o in ("-v", "--verbose"):
	    myopts.verbose = True
	    sys.stderr.write("Verbose output on\n")
    if len(args) < 1:
	usage()
	sys.exit(1)
    if myopts.verbose:
	sys.stderr.write("output = %s\n" % myopts.output)
	sys.stderr.write("dep-re = %s\n" % myopts.dep_re)
	sys.stderr.write("dep-template = %s\n" % myopts.dep_template)
	sys.stderr.write("mod-re = %s\n" % myopts.mod_re)
	sys.stderr.write("mod-template = %s\n" % myopts.mod_template)
	sys.stderr.write("o-prefix = %s\n" % myopts.o_prefix)
    return (myopts, args)

def main():
    (opts, filenames) = process_args()
    if opts.output:
	outf = open(opts.output, "w")
	if opts.verbose:
	    sys.stderr.write("Output to %s\n" % opts.output)
    else:
	outf = sys.stdout
	if opts.verbose:
	    sys.stderr.write("Output to STDOUT\n")
    outf.write("# DO NOT EDIT --- auto-generated file\n")

    modlist = {}

    for filename in filenames:
	(deps, mods) = get_deps_and_mods(filename, opts)

        filebase, fileext = os.path.splitext(filename)
        for cmod in mods:
            modlist[cmod] = opts.o_prefix+os.path.basename(filebase+".o")

    if opts.verbose:
        sys.stderr.write("modlist: %s\n" % " ".join(modlist))

    for filename in filenames:
	if opts.verbose:
	    sys.stderr.write("========= filename: %s\n" % filename)
	
	(deps, mods) = get_deps_and_mods(filename, opts)
	if deps:
            new_deps = []
            filebase, fileext = os.path.splitext(filename)
            objectfile = opts.o_prefix+os.path.basename(filebase+".o")
            for cdeps in deps:
                if modlist.get(cdeps, opts.o_prefix+cdeps) != objectfile:
                    if cdeps in modlist:
                        new_deps.append(modlist.get(cdeps, opts.o_prefix+cdeps))
	    write_deps(outf, filename, new_deps, opts)

	if opts.verbose:
            sys.stderr.write("deps: %s\n" % " ".join(deps))
	    sys.stderr.write("mods: %s\n" % " ".join(mods))
	    sys.stderr.write("newdeps: %s\n" % " ".join(new_deps))


    outf.close()

if __name__ == "__main__":
    main()
