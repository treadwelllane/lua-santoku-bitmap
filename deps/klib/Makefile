results.mk:
	wget https://github.com/attractivechaos/klib/raw/refs/heads/master/khash.h
	wget https://github.com/RoaringBitmap/CRoaring/releases/download/v4.3.1/roaring.c
	wget https://github.com/RoaringBitmap/CRoaring/releases/download/v4.3.1/roaring.h
	echo "LIB_CFLAGS += -I$(CURDIR)" >> results.mk
	touch $@
