package main

import (
	"math"
)


// assigns all atoms near heteroatoms to functional groups
func createFunctionalGroups(atoms map[int]*atom) int {

	hetAtomCounter := 0
	// assign innate heteroatoms
	for atomID, atom := range atoms {
		hetAtom := false
		// first all non H/C atoms are heteroatoms
		if atom.element != "C" && atom.element != "H" {
			hetAtom = true

		} else if atom.element == "C" && len(atom.bondedAtoms) < 4 {
			hetAtom = true
		}

		// if one is true, assign all atoms bonded to it to its group
		if hetAtom == true {
			hetAtomCounter++
			atom.isInFuncGroup = true
			for _, bondedAtom := range atom.bondedAtoms {
				if atoms[bondedAtom].element != "H" {
					atoms[bondedAtom].isInFuncGroup = true
					union(atoms, atomID, bondedAtom)
				}
			}
		}
	}

	return hetAtomCounter

}

// takes all adjacent atoms not in functional groups and merges them into groups
func mergeAlkanes(atoms map[int]*atom) {
	// Add all neighbors next to a heteroatom into its group
	//fmt.Println("Merging heteroatoms into functional groups...")
	for atomID, atom := range atoms {
		if atom.isInFuncGroup == false && atom.element == "C" {
			for _, bondedAtom := range atom.bondedAtoms {
				if atoms[bondedAtom].isInFuncGroup == false && atom.element == "C" {
					union(atoms, atomID, bondedAtom)
				}
			}
		}
	}
}

func mergeHydrogens(atoms map[int]*atom) {
	// add all hydrogens attached to functional group atoms to their functional group
	for atomID, atom := range atoms {
		if atom.element == "H" {
			union(atoms,atomID,atom.bondedAtoms[0])
		}
	}
}

// cleaves the bonds connecting groups
func separateGroups(atoms map[int]*atom) {
	for atomNum, atomInfo := range atoms {
		for _, bondedAtomNum := range atomInfo.bondedAtoms {
			if !connected(atoms,atomNum,bondedAtomNum) {
				removeBond(atoms,atomNum,bondedAtomNum)
			}
		}
	}
}

// removes bond between two atoms and adds hydrogens to molecule ends in its place
func removeBond(atoms map[int]*atom, atom1 int, atom2 int) {
	disconnect(atoms,atom1,atom2)
	addHydrogen(atoms,atom1,atom2)
	addHydrogen(atoms,atom2,atom1)
}

// removes each atom from the other's list of bonded atoms
func disconnect(atoms map[int]*atom, atom1 int, atom2 int) {
	for i, bondedAtom := range atoms[atom1].bondedAtoms {
		if bondedAtom == atom2 {
			atoms[atom1].bondedAtoms = remove(atoms[atom1].bondedAtoms,i)
		}
	}
	for i, bondedAtom := range atoms[atom2].bondedAtoms {
		if bondedAtom == atom1 {
			atoms[atom2].bondedAtoms = remove(atoms[atom2].bondedAtoms,i)
		}
	}
}

// removes element from slice
func remove(s []int, i int) []int {
	s[i] = s[len(s)-1]
	return s[:len(s)-1]
}

func addHydrogen(atoms map[int]*atom, atom1 int, atom2 int) {
	// calculate location of new h atom
	hpos := getHydrogenCoordinates(atoms,atom1,atom2)

	// create new h atom and set parameters
	var newAtom atom
	index := len(atoms)+1
	newAtom.element = "H"
	newAtom.pos = hpos
	newAtom.parent = index
	newAtom.atomType = 5
	newAtom.bondedAtoms = []int{atom1}
	newAtom.treeSize = 1
	newAtom.isInFuncGroup = true

	// adjust tree size of root
	atoms[root(atoms,atom1)].treeSize++

	// add to map

	atoms[index] = &newAtom
	// bind to atom1
	atoms[atom1].bondedAtoms = append(atoms[atom1].bondedAtoms, index)
	// join groups
	union(atoms,atom1,index)
}

func getHydrogenCoordinates(atoms map[int]*atom, atom1 int, atom2 int) []float64 {
	atom1coords := atoms[atom1].pos
	atom2coords := atoms[atom2].pos

	distance := 0.0
	for i := 0; i < 3; i++ {
		distance += math.Pow(atom2coords[i] - atom1coords[i],2)
	}
	distance = math.Sqrt(distance)

	r := make([]float64,3)
	for i:=0; i < 3; i++ {
		r[i] = (hydrogenBondDistance * (atom2coords[i] - atom1coords[i]) / distance) + atom1coords[i]
	}

	return r
}

// Creates a map of integer slices, where the key is the group identifier and the value contains
// the integer number of all atoms in a group together
func getGroups(atoms map[int]*atom) map[int][]int {
	// get groups
	groups := make(map[int][]int)
	for atomNum := range atoms {
		atomRoot := root(atoms,atomNum)
		// if group not in map already
		if _, ok := groups[atomRoot]; !ok {
			// add group to map
			groups[atomRoot] = []int{atomNum}
		} else {
			// append atom to group in map
			groups[atomRoot] = append(groups[atomRoot], atomNum)
		}
	}
	return groups
}

// groupMap connects the atom numbering within the old molecule (e.g. 1-100) to the atom numbering within the fragment (e.g. 1-10)
// groupSlice stores the
func groupMap(atoms map[int]*atom, group []int) (map[int]int, []int) {
	groupMap := make(map[int]int)
	groupSlice := make([]int,len(group))
	for i, atomID := range group {
		groupMap[atomID] = i+1
		groupSlice[i] = atomID
	}

	return groupMap, groupSlice
}


