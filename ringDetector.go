package main

func detectRings(atoms map[int]*atom) {
	// Copy from the original map to a temporary map
	atoms2 := make(map[int]*atom)
	for key, value := range atoms {
		atoms2[key] = value
	}
	// remove terminals
	removeTerminals(atoms2)
}

func removeTerminals(atoms2 map[int]*atom) {
	// track current terminals found this run through and last runthrough
	terminalsFoundLastIteration := -1
	terminalsFoundThisIteration := -1
	// need both to be zero as deleting terminals can create new terminals
	for terminalsFoundLastIteration == 0 && terminalsFoundThisIteration == 0 {
		terminalsFoundThisIteration = 0
		for atomID, atomInfo := range atoms2 {
			// find a terminal
			if len(atomInfo.bondedAtoms) < 2 {
				//
				terminalsFoundThisIteration++
				// remove references to this terminal in bonded atoms
				for i, neighborID := range atomInfo.bondedAtoms {
					atoms2[neighborID].bondedAtoms = remove(atoms2[neighborID].bondedAtoms, atomID)
				}
				// delete terminal from map
				delete(atoms2, atomID)

			}
		}

		terminalsFoundLastIteration = terminalsFoundThisIteration
	}
}

