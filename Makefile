all: clean addon

addon:
	zip -r qcbymm.zip qcbymm.r tool-inf LICENSE

clean:
	rm qcbymm.zip
