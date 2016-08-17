ACLOCAL_AMFLAGS = -I m4
SUBDIRS = frepattern functions lcm53 flowers
bin_PROGRAMS = lamp eliminate_comb flower
lamp_SOURCES = lamp.cpp LampCore.cpp LampCore.h ReadFile.cpp ReadFile.h Transaction.cpp Transaction.h
lamp_LDADD = frepattern/libfrepattern.a functions/libfunctions.a lcm53/liblcm.a \
             $(BOOST_FILESYSTEM_LIBS) $(BOOST_IOSTREAMS_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_SYSTEM_LIBS)
lamp_CPPFLAGS = $(BOOST_CPPFLAGS)
lamp_LDFLAGS = $(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_IOSTREAMS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_SYSTEM_LDFLAGS)
eliminate_comb_SOURCES = eliminate_comb.cpp
flower_SOURCES = flower.cpp LampCore.cpp LampCore.h ReadFile.cpp ReadFile.h Transaction.cpp Transaction.h
flower_LDADD = frepattern/libfrepattern.a functions/libfunctions.a lcm53/liblcm.a flowers/libflower.a \
             $(BOOST_FILESYSTEM_LIBS) $(BOOST_IOSTREAMS_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_SYSTEM_LIBS)
flower_CPPFLAGS = $(BOOST_CPPFLAGS)
flower_LDFLAGS = $(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_IOSTREAMS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_SYSTEM_LDFLAGS)

TESTS = test_lamp
check_PROGRAMS = test_lamp
test_lamp_SOURCES = test_lamp.cpp true_comb_list.h \
               LampCore.cpp LampCore.h ReadFile.cpp ReadFile.h Transaction.cpp Transaction.h
test_lamp_LDADD = frepattern/libfrepattern.a functions/libfunctions.a lcm53/liblcm.a \
             $(BOOST_FILESYSTEM_LIBS) $(BOOST_IOSTREAMS_LIBS) $(BOOST_SYSTEM_LIBS) $(GTEST_LIBS)
test_lamp_CPPFLAGS = $(BOOST_CPPFLAGS)
test_lamp_LDFLAGS = $(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_IOSTREAMS_LDFLAGS) $(BOOST_SYSTEM_LDFLAGS) $(GTEST_LDFLAGS)

EXTRA_DIST = sample

