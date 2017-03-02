'''
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
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
