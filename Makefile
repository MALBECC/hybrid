make:
	(cd src ; $(MAKE) )
	(cd tools ; $(MAKE) )
clean: 
	@echo "==> Cleaning object, library, and executable files"
	(cd src ; $(MAKE) clean)
	
