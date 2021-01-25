#########################
#
#    Bubblyhelix - clean.sh
#
#    Copyright 2012 Eivind Tostesen.
#
#    This file is part of Bubblyhelix.
#
#    Bubblyhelix is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Bubblyhelix is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Bubblyhelix.  If not, see <http://www.gnu.org/licenses/>.
#
#########################


# Removes all files beginning with .DS or ._ in the current directory and below
# Practical before creating tarball.

# usage: ./tools/clean.sh

# with prompting (y/n):
find . -name '.DS*' -ok rm {} \;
find . -name '._*' -ok rm {} \;

# without prompting:
#find . -name '.DS*' -exec rm {} \;
#find . -name '._*' -exec rm {} \;

