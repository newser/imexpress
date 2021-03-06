# Copyright (C) 2017 haniu (niuhao.cn@gmail.com)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
# USA.

# set an option's value and then hide it
macro(hide name type value)
    set(${name} ${value} CACHE ${type} "hidden" FORCE)
    mark_as_advanced(${name})
endmacro(hide)

# help to check if file exists
macro(check_exist result dir type name)
    set(${result} "")
    if (DEVELOP_BUILD)
        foreach(t ${type})
            foreach(n ${name})
                if (EXISTS ${dir}/${t}/${n})
                    set(${result} ${dir}/${t}/${n})
                    break()
                endif ()
            endforeach ()
        endforeach ()
    endif ()
endmacro(check_exist)
