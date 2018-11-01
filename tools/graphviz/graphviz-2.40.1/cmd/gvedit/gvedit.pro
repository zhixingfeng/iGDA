
LIBS += \
        -L$(top_builddir)/lib/gvc/.libs -lgvc \
        -L$(top_builddir)/lib/cgraph/.libs -lcgraph \
        -L$(top_builddir)/lib/cdt/.libs -lcdt \
        -lexpat -lz -lltdl

INCLUDEPATH += \
	../../lib/gvc \
	../../lib/common \
	../../lib/pathplan \
	../../lib/cgraph \
	../../lib/cdt \
	../..

CONFIG += qt
HEADERS = ../../cmd/gvedit/mainwindow.h ../../cmd/gvedit/mdichild.h ../../cmd/gvedit/csettings.h ../../cmd/gvedit/imageviewer.h ../../cmd/gvedit/ui_settings.h
SOURCES = ../../cmd/gvedit/main.cpp ../../cmd/gvedit/mainwindow.cpp ../../cmd/gvedit/mdichild.cpp ../../cmd/gvedit/csettings.cpp ../../cmd/gvedit/imageviewer.cpp
RESOURCES     = ../../cmd/gvedit/mdi.qrc

