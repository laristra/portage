'''
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
'''
import subprocess, sys
# This script parses the README file which describes how to build portage
# It looks for code blocks enclosed in ``` and requires that the first
# line of a block be a followed by a hostname.  The script will then
# execute all commands specified after that point on the host machine
# until the block is closed with a ```.
# Args: Arg 1 is file to parse, Arg 2 is the workspace
code = False
scripts = {}
mach_kw = "machine="
for line in open(sys.argv[1]):
    if line.startswith("```"):
        code = not code
        currentScript = None
    elif code:
        if currentScript is None:
            # we can have code, and the code can start with comments,
            # but we identify parts that should run on a particular machine
            # with the "machine=" keyword
            if not line.startswith("#") and mach_kw not in line: continue
            host = line.split(mach_kw)[-1].strip()
            currentScript = ["set -e\n pwd \n cd "+ sys.argv[2] + "\n rm -Rf build \n"]
            scripts[host] = currentScript
        else: currentScript.append(line)

retcode = 0
for host in scripts:
    if(host.startswith('varan')):
        print "HOST " + host
        scripts[host] += 'cd ..;rm -Rf build;'
        print scripts[host]
        sshHost = host.split(':')[0]
        proc = subprocess.Popen(['bash'], stdin = subprocess.PIPE)
        print "opened process"
        proc.communicate("".join(scripts[host]))
        if proc.returncode:
            retcode = proc.returncode
sys.exit(retcode)
