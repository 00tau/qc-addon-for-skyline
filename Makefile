all: clean addon

addon:
	zip -r easyQC.zip easyQC.r tool-inf LICENSE

clean:
	rm -f easyQC.zip
