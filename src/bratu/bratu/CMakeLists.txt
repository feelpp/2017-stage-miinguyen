###  CMakeLists.txt; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 29 Mar 2017
#
#  Copyright (C) 2017 Feel++ Consortium
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#


#feelpp_add_application( bratu_2d SRCS qs_bratu.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2)

#feelpp_add_application(bratu_2d SRCS bratu.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2)

#feelpp_add_application(bratu_3d SRCS bratu.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3)

#feelpp_add_application(bratu_mixte SRCS bratu_mixte.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2)

feelpp_add_application(bratu_fort_3d SRCS bratu_fort.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3)

#feelpp_add_application(bratu_fort_3D SRCS bratu_fort_3D.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3)

#feelpp_add_application(bratu_p_laplacian SRCS bratu_p_laplacian.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3)

#feelpp_add_application(bratu_p_laplacian_continu_3d SRCS bratu_p_laplacian_continu.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3)

#feelpp_add_application(bratu_p_laplacian_3d SRCS bratu_p_laplacian_3D.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3)


#feelpp_add_application(doc_nonlinearpow SRCS nonlinearpow.cpp LABELS doc)

#feelpp_add_application(doc_nlda SRCS nlda.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2)

#feelpp_add_application(doc_steadynavierstokes SRCS navierstokes.cpp CFG passerni.cfg kovasnay.cfg GEO kovasnaystraight2.geo 2D_FDA_L.geo LABELS doc )

#feelpp_add_application(doc_steadynavierstokes SRCS navierstokes.cpp CFG passerni.cfg GEO 2D_FDA_L.geo LABELS doc )
