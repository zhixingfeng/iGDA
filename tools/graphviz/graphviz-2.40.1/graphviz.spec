# $Id$ $Revision$
# graphviz.spec.  Generated from graphviz.spec.in by configure.

# Note: pre gd-2.0.34 graphviz uses its own gd tree with gif support and other fixes

#-- Global graphviz rpm and src.rpm tags-------------------------------------
Name:    graphviz
Summary: Graph Visualization Tools
Version: 2.40.1

%define truerelease 1
%{?distroagnostic: %define release %{truerelease}}
%{!?distroagnostic: %define release %{truerelease}%{?dist}}

Release: %{?release}
Group:   Applications/Multimedia
License: EPL
URL:     http://www.graphviz.org/
Source0: http://www.graphviz.org/pub/graphviz/development/SOURCES/graphviz-2.40.1.tar.gz

# graphviz is relocatable - Caution: this feature is used in AT&T,
#   but probably will not be supported in Redhat/Fedora/Centos distros
Prefix: /usr

#-- feature and package selection -------------------------------------------
#   depends on macro values set by redhat-rpm-config

# All features are off (undefined) by default
# To enable, use: <percent>define FEATURE 1
# Available features are:
#    SHARP GHOSTSCRIPT _GO GUILE _IO JAVA LUA OCAML ORTHO PERL PHP
#    PYTHON RUBY R_LANG TCL IPSEPCOLA MYLIBGD PANGOCAIRO RSVG
#    GTK GLITZ SMYRNA DEVIL MING GDK _QT WEBP

# SuSE uses a different mechanism to generate BuildRequires
# norootforbuild
# neededforbuild  expat freetype2-devel gcc tcl tcl-devel tk tk-devel x-devel-packages

BuildRoot:     %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)

BuildRequires: zlib-devel expat-devel ann-devel
BuildRequires: ksh bison m4 flex swig tk tcl >= 8.3 freetype-devel >= 2

#-- All platforms
%define __X 1

#-- Red Hat Enterprise Linux (also Centos) specific Build Requirements
%if 0%{?rhel}
%if %rhel == 4
BuildRequires: xorg-x11-devel
%endif
%if %rhel >= 4
BuildRequires: fontconfig-devel
%define TCL 1
BuildRequires: tcl-devel tk-devel 
%define PERL 1
%define RUBY 1
BuildRequires: ruby ruby-devel
%define GUILE 1
BuildRequires: guile-devel
%define PYTHON 1
BuildRequires: python python-devel
%define IPSEPCOLA 1
%define ORTHO 1
%endif
%if %rhel >= 5
BuildRequires: libXaw-devel libSM-devel libICE-devel libXpm-devel libXt-devel libXmu-devel libXext-devel libX11-devel
BuildRequires: libtool-ltdl libtool-ltdl-devel
%define JAVA 1
BuildRequires: java-devel
%define PHP 1
BuildRequires: php-devel
%define PANGOCAIRO 1
BuildRequires: cairo-devel >= 1.1.10 pango-devel gmp-devel gtk2-devel
%define RSVG 1
BuildRequires: librsvg2-devel
%define SFDP 1
BuildRequires: gts-devel
%define SMYRNA 1
BuildRequires: libglade2-devel gtkglarea2-devel gtkglext-devel freeglut-devel
%define GTK 1
%endif
%if %rhel >= 6
BuildRequires: perl-devel perl-libs perl-ExtUtils-Embed
%define LUA 1
BuildRequires: lua-devel
%define OCAML 1
BuildRequires: ocaml
%define MYLIBGD 0
BuildRequires: gd gd-devel
%define GDK 1
%define GHOSTSCRIPT 1
BuildRequires: ghostscript-devel
%define POPPLER 1
BuildRequires: poppler-glib-devel
%define _QT 1
BuildRequires: qt-devel
%endif
%endif

#-- Fedora specific Build Requirements
%if 0%{?fedora}
%if %fedora >= 9
BuildRequires: libXaw-devel libSM-devel libICE-devel libXpm-devel libXt-devel libXmu-devel libXext-devel libX11-devel
BuildRequires: fontconfig-devel
BuildRequires: libtool-ltdl libtool-ltdl-devel
%define TCL 1
BuildRequires: tcl-devel tk-devel 
%define PERL 1
BuildRequires: perl-devel perl-ExtUtils-Embed
%define RUBY 1
BuildRequires: ruby ruby-devel
%define GUILE 1
BuildRequires: guile-devel
%define PYTHON 1
BuildRequires: python python-devel
%define JAVA 1
BuildRequires: java-devel
%define _QT 1
BuildRequires: qt-devel
%ifnarch ppc64
%define SHARP 1
BuildRequires: mono-core
%define OCAML 1
BuildRequires: ocaml
%endif
%define LUA 1
BuildRequires: lua-devel
%define PANGOCAIRO 1
BuildRequires: cairo-devel >= 1.1.10 pango-devel gmp-devel gtk2-devel
%define RSVG 1
BuildRequires: librsvg2-devel
%define MYLIBGD 0
BuildRequires: gd gd-devel
%define DEVIL 1
BuildRequires: DevIL-devel
%define SFDP 1
BuildRequires: gts-devel
%define LASI 1
BuildRequires: lasi-devel
%define GDK 1
%define IPSEPCOLA 1
%define ORTHO 1
%define GTK 1
%define SMYRNA 1
BuildRequires: libglade2-devel gtkglarea2-devel gtkglext-devel glade3-libgladeui-devel freeglut-devel
#define GLITZ 1
# BuildRequires: glitz-devel
%define GHOSTSCRIPT 1
BuildRequires: ghostscript-devel
%define POPPLER 1
BuildRequires: poppler-glib-devel
#define MING 1
#BuildRequires: ming ming-devel
%define R_LANG 1
BuildRequires: R-devel swig >= 1.3.33 perl-libs
%endif
%if %fedora >= 16
%define WEBP 1
BuildRequires: libwebp-devel
%endif
%if %fedora >= 19
%define SEPARATE_LICENSE 1
%endif
%if %fedora >= 21
%define _GO 1
BuildRequires: golang >= 1.2 gcc-go swig >= 3.0.2
%endif
%if %fedora <= 24
%define PHP 1
BuildRequires: php-devel
# disable PHP7 in fc25 until swig supports it - requires swig-3.0.11
%endif
%endif

# What a meal PHP makes of versioning !!!
%define php_extdir %(php-config --extension-dir 2>/dev/null || echo %{_libdir}/php4)
%global php_apiver %((echo 0; php -i 2>/dev/null | sed -n 's/^PHP API => //p') | tail -1)

# Fix private-shared-object-provides
# RPM 4.8
%{?filter_provides_in: %filter_provides_in %{php_extdir}/.*\.so$}
%{?filter_setup}
# RPM 4.9
%global __provides_exclude_from %{?__provides_exclude_from:%__provides_exclude_from|}%{php_extdir}/.*\\.so$


#-- main graphviz rpm ------------------------------------------------
Requires:         %{name}-nox = %{version}-%{release}
Requires:         %{name}-x = %{version}-%{release}

%description
A collection of tools for the manipulation and layout
of graphs (as in nodes and edges, not as in bar-charts).

%files
%defattr(-,root,root,-)
%if 0%{?SEPARATE_LICENSE}
%license COPYING
%endif
%doc COPYING AUTHORS ChangeLog NEWS README

#-- graphviz-nox rpm --------------------------------------------------
%package nox
Group:            Applications/Multimedia
Summary:          Graphviz commands that do not depend on X11 - all graphviz installation will want this.
Requires:         urw-fonts
Requires:         %{name}-libs = %{version}-%{release}
Requires:         %{name}-plugins-core = %{version}-%{release}

%description nox
Graphviz plugins and commands that do not depend on x11. 
When used alone, this is intended for mimimalist web-server apps with no X11 on the server.

%files nox
%defattr(-,root,root,-)
%if 0%{?SEPARATE_LICENSE}
%license COPYING
%endif
%doc COPYING AUTHORS ChangeLog NEWS README
%exclude %{_bindir}/dot_builtins
%{_bindir}/acyclic
%{_bindir}/bcomps
%{_bindir}/ccomps
%{_bindir}/circo
%{_bindir}/cluster
%{_bindir}/dijkstra
%{_bindir}/dot
%{_bindir}/dot2gxl
%{_bindir}/fdp
%{_bindir}/gc
%{_bindir}/gml2gv
%{_bindir}/graphml2gv
%{_bindir}/gv2gml
%{_bindir}/gv2gxl
%{_bindir}/gvcolor
%{_bindir}/gvgen
%{_bindir}/gvmap
%{_bindir}/gvmap.sh
%{_bindir}/gvpack
%{_bindir}/gvpr
%{_bindir}/gxl2dot
%{_bindir}/gxl2gv
%{_bindir}/mingle
%{_bindir}/edgepaint
%{_bindir}/mm2gv
%{_bindir}/neato
%{_bindir}/nop
%{_bindir}/osage
%{_bindir}/patchwork
%{_bindir}/prune
%{_bindir}/sccmap
%{_bindir}/sfdp
%{_bindir}/tred
%{_bindir}/twopi
%{_bindir}/unflatten
%{_mandir}/man1/acyclic.1*
%{_mandir}/man1/bcomps.1*
%{_mandir}/man1/ccomps.1*
%{_mandir}/man1/circo.1*
%{_mandir}/man1/cluster.1*
%{_mandir}/man1/diffimg.1*
%{_mandir}/man1/dijkstra.1*
%{_mandir}/man1/dot.1*
%{_mandir}/man1/fdp.1*
%{_mandir}/man1/gc.1*
%{_mandir}/man1/gml2gv.1*
%{_mandir}/man1/graphml2gv.1*
%{_mandir}/man1/gv2gml.1*
%{_mandir}/man1/gv2gxl.1*
%{_mandir}/man1/gvcolor.1*
%{_mandir}/man1/gvgen.1*
%{_mandir}/man1/gvmap.1*
%{_mandir}/man1/gvmap.sh.1*
%{_mandir}/man1/gvpack.1*
%{_mandir}/man1/gvpr.1*
%{_mandir}/man1/gxl2gv.1*
%{_mandir}/man1/mingle.1*
%{_mandir}/man1/edgepaint.1*
%{_mandir}/man1/mm2gv.1*
%{_mandir}/man1/neato.1*
%{_mandir}/man1/nop.1*
%{_mandir}/man1/osage.1*
%{_mandir}/man1/patchwork.1*
%{_mandir}/man1/prune.1*
%{_mandir}/man1/sccmap.1*
%{_mandir}/man1/sfdp.1*
%{_mandir}/man1/tred.1*
%{_mandir}/man1/twopi.1*
%{_mandir}/man1/unflatten.1*
%{_mandir}/man7/graphviz.7*
%dir %{_datadir}/graphviz
%{_datadir}/graphviz/gvpr/*

#-- graphviz-libs rpm --------------------------------------------------
%package libs
Group:            Applications/Multimedia
Summary:          Graphviz base libs
Requires(post):   /sbin/ldconfig
Requires(postun): /sbin/ldconfig

%description libs
Graphviz core libs 

%files libs
%defattr(-,root,root,-)
%{_libdir}/libcdt.so.*
%{_libdir}/libcgraph.so.*
%{_libdir}/libgvc.so.*
%{_libdir}/libgvpr.so.*
%{_libdir}/libpathplan.so.*
%{_libdir}/libxdot.so.4*
%{_libdir}/liblab_gamut.so.*

#-- graphviz-plugins-core rpm --------------------------------------------------
%package plugins-core
Group:            Applications/Multimedia
Summary:          Graphviz plugins - core layout engines and text renderers
Requires:         %{name} = %{version}-%{release}

%description plugins-core
Graphviz plugins - core layout engines and text renderers

# run "dot -c" to generate plugin config in {_libdir}/graphviz/config
%post
LD_LIBRARY_PATH=$RPM_INSTALL_PREFIX0/%{_lib} $RPM_INSTALL_PREFIX0/bin/dot -c

# if there is no dot after everything else is done, then remove config
%postun
if [ $1 -eq 0 ]; then
        rm -f $RPM_INSTALL_PREFIX0/%{_lib}/graphviz/config || :
fi

%files plugins-core
%defattr(-,root,root,-)
%dir %{_libdir}/graphviz
%{_libdir}/graphviz/libgvplugin_core.so.*
%{_libdir}/graphviz/libgvplugin_dot_layout.so.*
%{_libdir}/graphviz/libgvplugin_neato_layout.so.*

#-- graphviz-x rpm --------------------------------------------------
%if 0%{?__X}
%package x
Group:            Applications/Multimedia
Summary:          Graphviz commands that depend on x11 - most installations will want this
Requires:         %{name}-nox = %{version}-%{release}
Requires:         %{name}-plugins-x = %{version}-%{release}

%description x
Graphviz plugins and commands that depend on x11 - most installations will want this

%files x
%defattr(-,root,root,-)
%if 0%{?SEPARATE_LICENSE}
%license COPYING
%endif
%doc COPYING AUTHORS ChangeLog NEWS README
%{_bindir}/lefty
%{_bindir}/lneato
%{_bindir}/dotty
%{_bindir}/vimdot
%{_mandir}/man1/lefty.1*
%{_mandir}/man1/lneato.1*
%{_mandir}/man1/dotty.1*
%{_mandir}/man1/vimdot.1*
%{_datadir}/graphviz/lefty
%if 0%{?SMYRNA}
%{_bindir}/smyrna
%{_datadir}/graphviz/smyrna
%{_mandir}/man1/smyrna.1*
%endif
%endif

#-- graphviz-plugins-x rpm --------------------------------------------------
%if 0%{?__X}
%package plugins-x
Group:            Applications/Multimedia
Summary:          Graphviz plugins and commands that depend on x11 - most installations will want this
Requires:         %{name} = %{version}-%{release}

%description plugins-x
Graphviz plugins that depend on x11 - most installations will want this.

# run "dot -c" to generate plugin config in {_libdir}/graphviz/config
%post plugins-x
LD_LIBRARY_PATH=$RPM_INSTALL_PREFIX0/%{_lib} $RPM_INSTALL_PREFIX0/bin/dot -c

%postun plugins-x
[ -x $RPM_INSTALL_PREFIX0/bin/dot ] && LD_LIBRARY_PATH=$RPM_INSTALL_PREFIX0/%{_lib} $RPM_INSTALL_PREFIX0/bin/dot -c || :

%files plugins-x
%defattr(-,root,root,-)
%dir %{_libdir}/graphviz
%{_libdir}/graphviz/libgvplugin_xlib.so.*
%if 0%{?PANGOCAIRO}
%{_libdir}/graphviz/libgvplugin_pango.so.*
%endif
%if 0%{?LASI}
%{_libdir}/graphviz/libgvplugin_lasi.so.*
%endif
%if 0%{?GDK}
%{_libdir}/graphviz/libgvplugin_gdk.so.*
%endif
%if 0%{?GTK}
%{_libdir}/graphviz/libgvplugin_gtk.so.*
%endif
%if 0%{?GHOSTSCRIPT}
%{_libdir}/graphviz/libgvplugin_gs.so.*
%endif
%if 0%{?POPPLER}
%{_libdir}/graphviz/libgvplugin_poppler.so.*
%endif
%if 0%{?RSVG}
%{_libdir}/graphviz/libgvplugin_rsvg.so.*
%endif
%endif

#-- graphviz-gd rpm --------------------------------------------------
%package gd
Group:            Applications/Multimedia
Summary:          Graphviz binaries that depend on libgd
Requires:         %{name} = %{version}-%{release}

# this next Requires is not strictly neccessary for diffimg,
#   but users will probably expect the plugins to get pulled in
Requires:         %{name}-plugins-gd = %{version}-%{release} 

%description gd
Graphviz binaries that depend on gd.  (Unless you absolutely have
to use GIF, you are recommended to use the PNG format instead because
of the better quality anti-aliased lines provided by the cairo+pango
based renderer.)

%files gd
%defattr(-,root,root,-)
%{_bindir}/diffimg

#-- graphviz-plugins-gd rpm --------------------------------------------------
%package plugins-gd
Group:            Applications/Multimedia
Summary:          Graphviz plugin for gd renderers.
Requires:         %{name} = %{version}-%{release}

%description plugins-gd
Graphviz plugin for image rendering using libgd.  (Unless you absolutely have
to use GIF, you are recommended to use the PNG format instead because
of the better quality anti-aliased lines provided by the cairo+pango
based renderer.)

# run "dot -c" to generate plugin config in {_libdir}/graphviz/config
%post plugins-gd
LD_LIBRARY_PATH=$RPM_INSTALL_PREFIX0/%{_lib} $RPM_INSTALL_PREFIX0/bin/dot -c

%postun plugins-gd
[ -x $RPM_INSTALL_PREFIX0/bin/dot ] && LD_LIBRARY_PATH=$RPM_INSTALL_PREFIX0/%{_lib} $RPM_INSTALL_PREFIX0/bin/dot -c || :

%files plugins-gd
%defattr(-,root,root,-)
%dir %{_libdir}/graphviz
%{_libdir}/graphviz/libgvplugin_gd.so.*

#-- graphviz-plugins-webp rpm --------------------------------------------------
%if 0%{?WEBP}
%package plugins-webp
Group:            Applications/Multimedia
Summary:          Graphviz plugin for webp format images, using libwebp
Requires:         %{name}-x = %{version}-%{release}
Obsoletes:        %{name}-webp

%description plugins-webp
Graphviz plugin for webp image rendering. 

# run "dot -c" to generate plugin config in {_libdir}/graphviz/config
%post plugins-webp
LD_LIBRARY_PATH=$RPM_INSTALL_PREFIX0/%{_lib} $RPM_INSTALL_PREFIX0/bin/dot -c

%postun plugins-webp
[ -x $RPM_INSTALL_PREFIX0/bin/dot ] && LD_LIBRARY_PATH=$RPM_INSTALL_PREFIX0/%{_lib} $RPM_INSTALL_PREFIX0/bin/dot -c || :

%files plugins-webp
%defattr(-,root,root,-)
%dir %{_libdir}/graphviz
%{_libdir}/graphviz/libgvplugin_webp.so.*
%endif

#-- graphviz-plugins-devil rpm --------------------------------------------------
%if 0%{?DEVIL}
%package plugins-devil
Group:            Applications/Multimedia
Summary:          Graphviz plugin for renderers based on DevIL
Requires:         %{name}-x = %{version}-%{release}
Obsoletes:        %{name}-devil

%description plugins-devil
Graphviz plugin for renderers based on DevIL.  (Unless you absolutely have
to use BMP, TIF, or TGA, you are recommended to use the PNG format instead
supported directly by the cairo+pango based renderer in the base graphviz rpm.)

# run "dot -c" to generate plugin config in {_libdir}/graphviz/config
%post plugins-devil
LD_LIBRARY_PATH=$RPM_INSTALL_PREFIX0/%{_lib} $RPM_INSTALL_PREFIX0/bin/dot -c

%postun plugins-devil
[ -x $RPM_INSTALL_PREFIX0/bin/dot ] && LD_LIBRARY_PATH=$RPM_INSTALL_PREFIX0/%{_lib} $RPM_INSTALL_PREFIX0/bin/dot -c || :

%files plugins-devil
%defattr(-,root,root,-)
%dir %{_libdir}/graphviz
%{_libdir}/graphviz/libgvplugin_devil.so.*
%endif

#-- graphviz-plugins-ming rpm --------------------------------------------------
%if 0%{?MING}
%package plugins-ming
Group:            Applications/Multimedia
Summary:          Graphviz plugin for flash renderer based on ming
Requires:         %{name}-x = %{version}-%{release}
Obsoletes:        %{name}-ming

%description plugins-ming
Graphviz plugin for -Tswf (flash) renderer based on ming.

# run "dot -c" to generate plugin config in {_libdir}/graphviz/config
%post plugins-ming
LD_LIBRARY_PATH=$RPM_INSTALL_PREFIX0/%{_lib} $RPM_INSTALL_PREFIX0/bin/dot -c

%postun plugins-ming
[ -x $RPM_INSTALL_PREFIX0/bin/dot ] && LD_LIBRARY_PATH=$RPM_INSTALL_PREFIX0/%{_lib} $RPM_INSTALL_PREFIX0/bin/dot -c || :

%files plugins-ming
%{_libdir}/graphviz/libgvplugin_ming.so.*
%{_libdir}/graphviz/*fdb
%endif

#-- graphviz-qt rpm --------------------------------------------------
%if 0%{?_QT}
%package qt
Group:            Applications/Multimedia
Summary:          Graphviz applications using _QT
Requires:         %{name}-x = %{version}-%{release}

%description qt
Graphviz applications using _QT - currently just gvedit

%files qt
%defattr(-,root,root,-)
%{_bindir}/gvedit
%{_datadir}/graphviz/gvedit
%{_mandir}/man1/gvedit.1*
%endif

#-- graphviz-lang-sharp rpm --------------------------------------------
%if 0%{?SHARP}
%package lang-sharp
Group:          Applications/Multimedia
Summary:        C# extension for graphviz
Requires:       %{name} = %{version}-%{release}, mono-core
Obsoletes:	%{name}-sharp

%description lang-sharp
C# extension for graphviz.

%files lang-sharp
%defattr(-,root,root,-)
%dir %{_libdir}/graphviz/sharp
%{_libdir}/graphviz/sharp/*
%{_mandir}/man3/*.3sharp.*
%endif

#-- graphviz-lang-go rpm --------------------------------------------
%if 0%{?_GO}
%package lang-go
Group:          Applications/Multimedia
Summary:        GO extension for graphviz
Requires:       %{name} = %{version}-%{release}, golang
Obsoletes:	%{name}-go

%description lang-go
Guile extension for graphviz.

%files lang-go
%defattr(-,root,root,-)
%dir %{_libdir}/graphviz/go
%{_libdir}/graphviz/go/*
%{_mandir}/man3/*.3go.*
%endif

#-- graphviz-lang-guile rpm --------------------------------------------
%if 0%{?GUILE}
%package lang-guile
Group:          Applications/Multimedia
Summary:        Guile extension for graphviz
Requires:       %{name} = %{version}-%{release}, guile
Obsoletes:	%{name}-guile

%description lang-guile
Guile extension for graphviz.

%files lang-guile
%defattr(-,root,root,-)
%dir %{_libdir}/graphviz/guile
%{_libdir}/graphviz/guile/*
%{_mandir}/man3/*.3guile.*
%endif

#-- graphviz-lang-io rpm -----------------------------------------------
%if 0%{?_IO}
%package lang-io
Group:          Applications/Multimedia
Summary:        Io extension for graphviz
Requires:       %{name} = %{version}-%{release}, io
Obsoletes:	%{name}-io

%description lang-io
Io extension for graphviz.

%files lang-io
%defattr(-,root,root,-)
%{_mandir}/man3/*.3io.*
%endif

#-- graphviz-lang-java rpm ---------------------------------------------
%if 0%{?JAVA}
%package lang-java
Group:          Applications/Multimedia
Summary:        Java extension for graphviz
Requires:       %{name} = %{version}-%{release}, java
Obsoletes:	%{name}-java

%description lang-java
Java extension for graphviz.

%files lang-java
%defattr(-,root,root,-)
%dir %{_libdir}/graphviz/java
%{_libdir}/graphviz/java/*
%{_mandir}/man3/*.3java.*
%endif

#-- graphviz-lang-lua rpm ----------------------------------------------
%if 0%{?LUA}
%package lang-lua
Group:          Applications/Multimedia
Summary:        Lua extension for graphviz
Requires:       %{name} = %{version}-%{release}, lua
Obsoletes:	%{name}-lua

%description lang-lua
Lua extension for graphviz.

%files lang-lua
%defattr(-,root,root,-)
%{_libdir}/lua*/*
%{_datadir}/graphviz/demo/*.lua*
%{_mandir}/man3/*.3lua.*
%exclude %{_libdir}/graphviz/lua/*.so
%endif

#-- graphviz-lang-ocaml rpm --------------------------------------------
%if 0%{?OCAML}
%package lang-ocaml
Group:          Applications/Multimedia
Summary:        Ocaml extension for graphviz
Requires:       %{name} = %{version}-%{release}, ocaml
Obsoletes:	%{name}-ocaml

%description lang-ocaml
Ocaml extension for graphviz.

%files lang-ocaml
%defattr(-,root,root,-)
%dir %{_libdir}/graphviz/ocaml
%{_libdir}/graphviz/ocaml/*
%{_mandir}/man3/*.3ocaml.*
%endif

#-- graphviz-lang-perl rpm ---------------------------------------------
%if 0%{?PERL}
%package lang-perl
Group:          Applications/Multimedia
Summary:        Perl extension for graphviz
Requires:       %{name} = %{version}-%{release}, perl
Obsoletes:	%{name}-perl

%description lang-perl
Perl extension for graphviz.

%files lang-perl
%defattr(-,root,root,-)
%{_libdir}/perl*/*
%{_datadir}/graphviz/demo/*.pl*
%{_mandir}/man3/*.3perl.*
%exclude %{_libdir}/graphviz/perl/*.so
%exclude %{_libdir}/graphviz/perl/*.pm
%endif

#-- graphviz-lang-php rpm ----------------------------------------------
%if 0%{?PHP}
%package lang-php
Group:          Applications/Multimedia
Summary:        PHP extension for graphviz
Requires:       %{name} = %{version}-%{release}, /usr/bin/php
Obsoletes:	%{name}-php

%description lang-php
PHP extension for graphviz.

%files lang-php
%defattr(-,root,root,-)
%config(noreplace) %{_sysconfdir}/php.d/%{name}.ini
%{php_extdir}/gv.so
%{_datadir}/php*/*
%{_datadir}/graphviz/demo/*.php*
%{_mandir}/man3/*.3php.*
%exclude %{_libdir}/graphviz/php/*.so
%exclude %{_libdir}/graphviz/php/*.php
%endif

#-- graphviz-lang-python rpm -------------------------------------------
%if 0%{?PYTHON}
%package lang-python
Group:          Applications/Multimedia
Summary:        Python extension for graphviz
Requires:       %{name} = %{version}-%{release}, python
Obsoletes:	%{name}-python

%description lang-python
Python extension for graphviz.

%files lang-python
%defattr(-,root,root,-)
%{_libdir}/python*/*
%{_datadir}/graphviz/demo/*.py*
%{_mandir}/man3/*.3python.*
%exclude %{_libdir}/graphviz/python/*.so
%exclude %{_libdir}/graphviz/python/*.py*
%endif

#-- graphviz-lang-R rpm ---------------------------------------------
%if 0%{?R_LANG}
%package lang-R
Group:          Applications/Multimedia
Summary:        R extension for graphviz
Requires:       %{name} = %{version}-%{release}, R
Obsoletes:	%{name}-R

%description lang-R
R extension for graphviz.

%files lang-R
%defattr(-,root,root,-)
%dir %{_libdir}/graphviz/R
%{_libdir}/graphviz/R/*
%{_mandir}/man3/*.3r.*
%endif

#-- graphviz-lang-ruby rpm ---------------------------------------------
%if 0%{?RUBY}
%package lang-ruby
Group:          Applications/Multimedia
Summary:        Ruby extension for graphviz
Requires:       %{name} = %{version}-%{release}, ruby
Obsoletes:	%{name}-ruby

%description lang-ruby
Ruby extension for graphviz.

%files lang-ruby
%defattr(-,root,root,-)
%{_libdir}/*ruby*/*
%{_datadir}/graphviz/demo/*.rb*
%{_mandir}/man3/*.3ruby.*
%exclude %{_libdir}/graphviz/ruby/*.so
%endif

#-- graphviz-lang-tcl rpm ----------------------------------------------
%if 0%{?TCL}
%package lang-tcl
Group:          Applications/Multimedia
Summary:        Tcl extension & tools for graphviz
Requires:       %{name} = %{version}-%{release}, tcl >= 8.3
Obsoletes:	%{name}-tcl

%description lang-tcl
Various tcl packages (extensions) for the graphviz tools.

%files lang-tcl
%defattr(-,root,root,-)
%{_libdir}/tcl*/*
%{_datadir}/graphviz/demo/*.tcl*
%{_datadir}/graphviz/demo/*_data
%{_datadir}/graphviz/demo/entities.html
%{_mandir}/man3/*.3tcl.*
%exclude %{_libdir}/graphviz/tcl/*
%endif

#-- graphviz-devel rpm --------------------------------------------
%package devel
Group:          Development/Libraries
Summary:        Development package for graphviz
Requires:       %{name}-libs = %{version}-%{release}, pkgconfig

%description devel
A collection of tools for the manipulation and layout
of graphs (as in nodes and edges, not as in bar-charts).
This package contains development files for graphviz-libs.

%files devel
%defattr(-,root,root,-)
%{_includedir}/graphviz
%{_libdir}/libcdt.so
%{_mandir}/man3/cdt.3.*
%{_libdir}/pkgconfig/libcdt.pc
%{_libdir}/libcgraph.so
%{_mandir}/man3/cgraph.3.*
%{_libdir}/pkgconfig/libcgraph.pc
%{_libdir}/libgvc.so
%{_mandir}/man3/gvc.3.*
%{_libdir}/pkgconfig/libgvc.pc
%{_libdir}/libgvpr.so
%{_mandir}/man3/gvpr.3.*
%{_libdir}/pkgconfig/libgvpr.pc
%{_libdir}/libpathplan.so
%{_mandir}/man3/pathplan.3.*
%{_libdir}/pkgconfig/libpathplan.pc
%{_libdir}/libxdot.so
%{_mandir}/man3/xdot.3.*
%{_libdir}/pkgconfig/libxdot.pc
%{_libdir}/liblab_gamut.so
%{_mandir}/man3/lab_gamut.3.*
%{_libdir}/pkgconfig/liblab_gamut.pc
%{_mandir}/man3/expr.3.*
%{_mandir}/man3/pack.3.*
%exclude %{_libdir}/graphviz/libgvplugin*
%exclude %{_libdir}/graphviz/*.so

#-- graphviz-graphs rpm -------------------------------------------
%package graphs
Group:          Applications/Multimedia
Summary:        Demo graphs for graphviz
%if 0%{?fedora} >= 11
BuildArch: noarch
%endif
%if 0%{?rhel} >= 6
BuildArch: noarch
%endif

%description graphs
Some demo graphs for graphviz.

%files graphs
%defattr(-,root,root,-)
%dir %{_datadir}/graphviz
%{_datadir}/graphviz/graphs
%if 0%{?SMYRNA}
%{_datadir}/graphviz/examples
%endif

#-- graphviz-doc rpm ----------------------------------------------
%package doc
Group:          Documentation
Summary:        PDF and HTML documents for graphviz
%if 0%{?fedora} >= 11
BuildArch: noarch
%endif
%if 0%{?rhel} >= 6
BuildArch: noarch
%endif

%description doc
Provides some additional PDF and HTML documentation for graphviz.

%files doc
%defattr(-,root,root,-)
%doc __doc/*

#-- building --------------------------------------------------

%prep
%setup -q
### the rpm in centos4/5 is too old to handle this
#for p in %{patches}; do
#     patch -p1 < ${p} || exit 1
#done

%build
# XXX ix86 only used to have -ffast-math, let's use everywhere
%{expand: %%define optflags %{optflags} -ffast-math}

# %%configure is broken in RH7.3 rpmbuild
CFLAGS="$RPM_OPT_FLAGS" \
./configure \
        --prefix=%{_prefix} \
        --bindir=%{_bindir} \
        --libdir=%{_libdir} \
        --includedir=%{_includedir} \
        --datadir=%{_datadir} \
        --mandir=%{_mandir} \
        --disable-static \
        --disable-dependency-tracking \
        --enable-sharp%{!?SHARP:=no} \
        --enable-go%{!?_GO:=no} \
        --enable-guile%{!?GUILE:=no} \
        --enable-io%{!?_IO:=no} \
        --enable-java%{!?JAVA:=no} \
        --enable-lua%{!?LUA:=no} \
        --enable-ocaml%{!?OCAML:=no} \
        --enable-perl%{!?PERL:=no} \
        --enable-php%{!?PHP:=no} \
        --enable-python%{!?PYTHON:=no} \
        --enable-r%{!?R_LANG:=no} \
        --enable-ruby%{!?RUBY:=no} \
        --enable-tcl%{!?TCL:=no} \
        --with%{!?DEVIL:out}-devil \
        --with%{!?WEBP:out}-webp \
        --with%{!?GDK:out}-gdk \
        --with%{!?GHOSTSCRIPT:out}-ghostscript \
        --with%{!?GLITZ:out}-glitz \
        --with%{!?IPSEPCOLA:out}-ipsepcola \
        --with%{!?LASI:out}-lasi \
        --with%{!?MING:out}-ming \
        --with%{!?_QT:out}-qt \
        --with%{!?PANGOCAIRO:out}-pangocairo \
        --with%{!?POPPLER:out}-poppler \
        --with%{!?RSVG:out}-rsvg \
        --with%{!?ORTHO:out}-ortho \
        --with%{!?SFDP:out}-sfdp \
        --with%{!?SMYRNA:out}-smyrna \
        --with%{!?__X:out}-x
make %{?_smp_mflags}

%install
rm -rf %{buildroot} __doc
make DESTDIR=%{buildroot} \
        docdir=%{buildroot}%{_docdir}/%{name} \
        pkgconfigdir=%{_libdir}/pkgconfig \
        install
find %{buildroot} -type f -name "*.la" -exec rm -f {} ';'
chmod -x %{buildroot}%{_datadir}/%{name}/lefty/*
cp -a %{buildroot}%{_datadir}/%{name}/doc __doc
rm -rf %{buildroot}%{_datadir}/%{name}/doc
%if 0%{?PHP}
# PHP configuration file
%{__mkdir_p} %{buildroot}%{_sysconfdir}/php.d
%{__cat} << __EOF__ > %{buildroot}%{_sysconfdir}/php.d/%{name}.ini
; Enable %{name} extension module
extension=gv.so
__EOF__
%endif

%check
%ifnarch ppc64 ppc
# regression test, segfaults on ppc/ppc64, possible endian issues?
# regressions tests require ksh93 - not available on centos3
#cd rtest
#make rtest
%endif

%clean
# optional regression test using the products in the build tree
%if 0%{?rtest}
cd rtest
make rtest
%endif
# clean up temporary installation
rm -rf %{buildroot}

#-- changelog --------------------------------------------------

%changelog
* Tue Mar 03 2009 John Ellson <ellson@research.att.com>
- release graphviz-2.22 - see ChangeLog for details

* Fri Oct 26 2007 John Ellson <ellson@research.att.com>
- reinstated binary relocatability in AT&T's rpms since we need it internally
* Wed Aug 15 2007 John Ellson <ellson@research.att.com>
- release 2.14.1 - see ChangeLog for details
* Thu Aug 2 2007 John Ellson <ellson@research.att.com>
- release 2.14 - see ChangeLog for details
* Fri Mar 16 2007 Stephen North <north@research.att.com>
- remove xorg-X11-devel from rhel >= 5
* Mon Dec 11 2006 John Ellson <john.ellson@comcast.net>
- fix graphviz-lua description (Fedora BZ#218191)
* Tue Sep 13 2005 John Ellson <ellson@research.att.com>
- split out language bindings into their own rpms so that 
  main rpm doesn't depend on (e.g.) ocaml

* Sat Aug 13 2005 John Ellson <ellson@research.att.com>
- imported various fixes from the Fedora-Extras .spec by Oliver Falk <oliver@linux-kernel.at>

* Wed Jul 20 2005 John Ellson <ellson@research.att.com>
- release 2.4
