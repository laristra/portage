'''
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
'''
import subprocess, sys
#This script parses the README file which describes how to build portage
#It looks for code blocks enclosed in ``` and requires that the first
#line of a block be a followed by a hostname.  The script will then
#execute all commands specified after that point on the host machine
#until the block is closed with a ```.
#Args: Arg 1 is file to parse, Arg 2 is the workspace
code = False
scripts = {}
for line in open(sys.argv[1]):
    if line.startswith("```"):
        code = not code
        currentScript = None
    elif code:
        if currentScript is None:
            assert line.startswith("#")
            currentScript = ["set -e\n pwd \n cd "+ sys.argv[2] + "\n rm -Rf build \n"]
            scripts[line.lstrip("#").strip()] = currentScript
        else: currentScript.append(line)

retcode = 0
for host in scripts:
    if(host.startswith('varan')):
        print "HOST " + host
        scripts[host] += 'cd ..;rm -Rf build;'
        print scripts[host]
        sshHost = host.split(':')[0]
        proc = subprocess.Popen(['ssh', '-T', "jenksrvngc@" + sshHost],
                                stdin = subprocess.PIPE)
        print "opened process"
        proc.communicate("".join(scripts[host]))
        if proc.returncode:
            retcode = proc.returncode
sys.exit(retcode)
