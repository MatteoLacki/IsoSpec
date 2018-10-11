* Properly name the classes:
	* conf -> subisotopologue
	* iso  -> molecule
	* partialExpProbs -> partialProbs
	* L... -> log_...

* Types:
	* why get_isotopeNo is an int and not a signed int?
	* the same with the subisotopologue configurations.


* general WTF:
	* get_conf_signature ??? 
	* IsoOrderedGenerator::candidate ???
	* IsoOrderedGenerator::terminate_search ???
	* IsoLayeredGenerator::advanceToNextConfiguration
		* why does it print out things to the console???
		* why is advanceToNextConfiguration_internal public, if it's called internal?
			* it ain't used anywhere else beyond that class anyway.

