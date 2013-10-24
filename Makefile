# makefile for plinktools by Iain Bancarz, ib5@sanger.ac.uk

PREFIX = 	PREFIX_DIRECTORY    # dummy value as default
DEST =		$(PREFIX)/plinktools

usage:
	@echo -e "Usage: make install PREFIX=<destination directory>\nWill install to the plinktools subdirectory of PREFIX.\nPREFIX must exist, plinktools subdirectory will be created if necessary."

install: $(PREFIX)
	install -d $(DEST)
	@rm -f $(DEST)/*.pyc    # force recompiling of any .pyc files
	install checksum.py compare.py het_by_maf.py __init__.py merge_bed.py plink.py $(DEST)
	@echo -e "Plinktools successfully installed."
