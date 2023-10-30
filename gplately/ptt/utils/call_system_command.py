
"""
    Copyright (C) 2017 The University of Sydney, Australia
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

from __future__ import print_function

import subprocess
import sys


# Function to call a command on the system (based on 'subprocess' module).
#   Option to check return code of command and which return code to check (defaults to checking code 0 for success).
#   Option to raise error or return False on failure (raises by default).
#   Option to print an error message to stderr on failure (prints by default).
#   Option to pass a stdin string to the standard input of command (not passed by default).
#   Option to receive stdout/stderr strings from the standard output/error of command (not received by default).
#   Optionally pass advanced parameters to 'subprocess.Popen()' such as 'shell=True' (no advanced parameters passed by default).
#
# On success:
#    Returns (stdout, stderr) tuple of strings if 'return_stdout' and 'return_stderr' are True,
#    else returns stdout string if only 'return_stdout' is True,
#    else returns stderr string if only 'return_stderr' is True,
#    else returns True (if both 'return_stdout' and 'return_stderr' are False).
#
# On failure:
#    Raises an exception if 'raise_errors' is True,
#    else returns None (which can be tested as if it was 'False').
#
def call_system_command(
        args,  # Command and its arguments - either a single string or a sequence of arguments (see subprocess.Popen()).
        check_return_code=0,  # Check command's return code with this value (set to None to avoid checking).
        raise_errors=True,  # Whether to raise an exception when there's an error (terminates calling script unless caught).
        print_errors=True,  # Whether to print an error message to stderr when there's an error.
        stdin=None,  # Optional string to send to stdin of the command.
        return_stdout=False,  # Whether to capture, and return, stdout of the command,
        return_stderr=False,  # Whether to capture, and return, stderr of the command,
        **subprocess_options):  # Advanced options passed directly to subprocess.Popen().
    
    # Whether to send stdin to command.
    if stdin is not None:
        stdin_pipe = subprocess.PIPE
    else:
        stdin_pipe = None
    
    # Whether to receive stdout/stderr from command.
    if return_stdout:
        stdout_pipe = subprocess.PIPE
    else:
        stdout_pipe = None
    if return_stderr:
        stderr_pipe = subprocess.PIPE
    else:
        stderr_pipe = None
    
    # Execute command.
    try:
        command = subprocess.Popen(args, stdin=stdin_pipe, stdout=stdout_pipe, stderr=stderr_pipe, universal_newlines=True, **subprocess_options)
        stdout, stderr = command.communicate(stdin)
    except ValueError as e:
        if print_errors:
            print("System command called with invalid arguments: {0}".format(e), file=sys.stderr)
        if not raise_errors:
            return None
        raise
    except OSError as e:
        if print_errors:
            print("Unable to execute system command: {0} {1}".format(args, e), file=sys.stderr)
        if not raise_errors:
            return None
        raise
    
    # Check return code (if requested).
    if check_return_code is not None:
        command_return_code = command.poll()
        if command_return_code != check_return_code:
            if print_errors:
                print("System command failed: {0} return code: {1}".format(args, command_return_code), file=sys.stderr)
            if not raise_errors:
                return None
            
            # Raise same error that subprocess.check_call() does.
            raise subprocess.CalledProcessError(command_return_code, args, output=stdout)
    
    # Return (stdout, stderr) tuple if requested.
    if return_stdout and return_stderr:
        return stdout, stderr
    elif return_stdout:
        return stdout
    elif return_stderr:
        return stderr
    
    return True


#   if __name__ == '__main__':
#
#       # Windows 'dir' command (change to 'ls -l', for example, on Mac and Linux).
#       call_system_command('dir /w', shell=True)
#       call_system_command(['dir', '/w'], print_errors=False, shell=True)
#       if not call_system_command(['dir', '/w'], raise_errors=False, shell=True):
#           print 'Failed directory listing.'
#
#       # Info on a grid file.
#       call_system_command(['gmt', 'grdinfo', 'mohoGOCE180.nc'])
#
#       # Sample a grid file at specified points.
#       # Send stdin (the lon/lat points to sample) and receive stdout (the sampled raster values).
#       output = call_system_command(["gmt", "grdtrack", "-nl", "-GmohoGOCE180.nc"], stdin='10 20\n-10 -20\n', return_stdout=True, return_stderr=True, raise_errors=False)
#       # Check return value since we're not raising exceptions (due to "raise_errors=False").
#       # Note: Need to test 'output' *before* attempting to extract tuple components from it.
#       if output:
#           stdout, stderr = output
#           print stdout
#       # Run again but only return stdout (not stderr).
#       stdout = call_system_command(["gmt", "grdtrack", "-nl", "-GmohoGOCE180.nc"], stdin='10 20\n-10 -20\n', return_stdout=True, raise_errors=False)
#       # Check return value since we're not raising exceptions (due to "raise_errors=False").
#       if stdout is not None:
#           print stdout
#           pass
