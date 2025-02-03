#
#    Copyright (C) 2024 The University of Sydney, Australia
#
#    This program is free software; you can redistribute it and/or modify it under
#    the terms of the GNU General Public License, version 2, as published by
#    the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

"""
Imports [`pyGPlates`](https://www.gplates.org/docs/pygplates/index.html).

Note that an alternative to `import gplates.pygplates` is just `import pygplates`.

This module was initially provided to support pickling of some pygplates classes.
But as of pygplates version 1.0, pickling is natively supported (within pygplates).
So now this module is no longer necessary, but retained for backward compatibility.
"""

from pygplates import *
