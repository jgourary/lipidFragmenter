package main

import (
	"errors"
	"log"
	"strconv"
)

///////////////////
// Atom data structure
///////////////////

type atom struct {
	element string
	atomType int
	bondedAtoms []int

	pos []float64

	// used for union-find algorithm
	parent int
	treeSize int
	isInFuncGroup bool

	// used for ring detection algorithm
	visited bool
	discTime int
	minDiscTime int
	parentBF int
	isCyclic bool

	// used for fragment to dimer recombination (this process also uses the ring detection vars)
	isTerminus bool
	isAnteTerminus bool
}

func copyAtom(thisAtom *atom) atom {
	var newAtom atom


	newAtom.element = thisAtom.element
	newAtom.atomType = thisAtom.atomType
	newAtom.bondedAtoms = make([]int, len(thisAtom.bondedAtoms))
	for i := range thisAtom.bondedAtoms {
		newAtom.bondedAtoms[i] = thisAtom.bondedAtoms[i]
	}

	newAtom.pos = thisAtom.pos

	// used for union-find algorithm
	newAtom.parent = thisAtom.parent
	newAtom.treeSize = thisAtom.treeSize
	newAtom.isInFuncGroup = thisAtom.isInFuncGroup

	// used for ring detection algorithm
	newAtom.visited = thisAtom.visited
	newAtom.discTime  = thisAtom.discTime
	newAtom.minDiscTime = thisAtom.discTime
	newAtom.parentBF = thisAtom.parentBF
	newAtom.isCyclic = thisAtom.isCyclic

	// used for fragment to dimer recombination (this process also uses the ring detection vars)
	newAtom.isTerminus = thisAtom.isTerminus
	newAtom.isAnteTerminus  = thisAtom.isAnteTerminus

	return newAtom
}

func copyMolecule(atomsA map[int]*atom) map[int]*atom {
	atomsB := make(map[int]*atom)
	for atomID, thisAtom := range atomsA {
		newAtom := copyAtom(thisAtom)
		atomsB[atomID] = &newAtom
	}
	return atomsB
}

// Includes renumbering atom map keys from 1 to len(molecule)-1
func copyMoleculeWithRenumbering(atomsA map[int]*atom) map[int]*atom {

	atomIDOldToNewMap := make(map[int]int)
	i := 1
	for oldAtomID := range atomsA {
		atomIDOldToNewMap[oldAtomID] = i
		i++
	}

	atomsB := make(map[int]*atom)
	for oldIndex, thisAtom := range atomsA {
		newAtom := copyAtom(thisAtom)
		newIndex := atomIDOldToNewMap[oldIndex]
		for j := 0; j < len(newAtom.bondedAtoms); j++ {
			newAtom.bondedAtoms[j] = atomIDOldToNewMap[newAtom.bondedAtoms[j]]
		}
		newAtom.parent = atomIDOldToNewMap[newAtom.parent]
		newAtom.parentBF = atomIDOldToNewMap[newAtom.parentBF]
		atomsB[newIndex] = &newAtom
	}
	return atomsB
}

// /////////////////
// Functions
// /////////////////

// assigns all atoms near heteroatoms to functional groups
func createFunctionalGroups(atoms map[int]*atom) int {

	hetAtomCounter := 0
	// assign innate heteroatoms
	for atomID, atom := range atoms {
		hetAtom := false
		// first, all non H/C atoms are heteroatoms
		if atom.element != "C" && atom.element != "H" {
			hetAtom = true
			// second, all C atoms with pi bonds are heteroatoms
		} else if atom.element == "C" && len(atom.bondedAtoms) < 4 {
			hetAtom = true
			// third, all atoms without a bridge type connection are part of a cycle and thus heteroatoms
		} else if atom.isCyclic == true {
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

// Creates a map of integer slices, where the key is the group identifier
func getRoots(atoms map[int]*atom) []int {
	// get roots
	roots := make(map[int]int)
	for atomNum := range atoms {
		atomRoot := root(atoms,atomNum)
		// if group not in map already
		if _, ok := roots[atomRoot]; !ok {
			// add group to map
			roots[atomRoot] = atomRoot
		}
	}
	rootsSlice := make([]int, len(roots))
	i := 0
	for k := range roots {
		rootsSlice[i] = k
		i++
	}
	//fmt.Println("debug: roots = ")
	//fmt.Println(rootsSlice)
	return rootsSlice
}

func getSingleFragments(atoms map[int]*atom, borderBonds [][2]int) ([]map[int]*atom, []bool) {

	roots := getRoots(atoms)
	singleFragSlice := make([]map[int]*atom, len(roots))
	isFragHydrocarbon := make([]bool, len(roots))

	// iterate through all border bonds
	for i := 0; i < len(borderBonds); i++ {
		// disconnect all bonds connecting fragments
		removeBondAndCapEnds(atoms,borderBonds[i][0],borderBonds[i][1])
	}

	// iterate through all roots
	for i := 0; i < len(roots); i++ {
		// create a copyFile of the molecule
		moleculeCopy := copyMolecule(atoms)

		// delete all atoms not of that root
		deleteAtomsNotOfGroup(moleculeCopy, []int {roots[i]})

		isHydrocarbon := true
		for atomID := range moleculeCopy {
			if moleculeCopy[atomID].element != "H" && moleculeCopy[atomID].element != "C" {
				isHydrocarbon = false
			}
		}
		//fmt.Println("debug: len of fragment = " + strconv.Itoa(len(moleculeCopy)))
		// copyFile remaining fragment with renumbering from 1...n
		singleFragSlice[i] = copyMoleculeWithRenumbering(moleculeCopy)
		isFragHydrocarbon[i] = isHydrocarbon
	}

	return singleFragSlice, isFragHydrocarbon
}

// returns all bonds connecting two different fragments
func getFragmentBorderBonds(atoms map[int]*atom) [][2]int {

	var borderBonds [][2]int

	for atomNum, atomInfo := range atoms {
		for _, bondedAtomNum := range atomInfo.bondedAtoms {
			if !connected(atoms,atomNum,bondedAtomNum) {
				// avoid counting the same bond twice
				if atomNum < bondedAtomNum {
					thisBond := [2]int {atomNum,bondedAtomNum}
					borderBonds = append(borderBonds,thisBond)
				}
			}
		}
	}
	return borderBonds
}



func getDoubleFragments(atoms map[int]*atom, borderBonds [][2]int) ([]map[int]*atom, []map[int]*atom, [][2]bool) {

	doubleFragSlice := make([]map[int]*atom, len(borderBonds))
	dimersSlice := make([]map[int]*atom, len(borderBonds))
	doubleFragIsHydrocarbon := make([][2]bool, len(borderBonds))

	// iterate through all border bonds
	for i := 0; i < len(borderBonds); i++ {
		// create copyFile of the molecule
		moleculeCopy := copyMolecule(atoms)

		// note the groups of the two fragments linked by the border bond
		group1 := root(moleculeCopy, borderBonds[i][0])
		group2 := root(moleculeCopy, borderBonds[i][1])


		//fmt.Println("debug: undoing border bonds and capping ends")
		// undo all border bonds except the one connecting the double fragment
		for j := 0; j < len(borderBonds); j++ {
			if i != j {
				removeBondAndCapEnds(moleculeCopy,borderBonds[j][0],borderBonds[j][1])
			}
		}

		//fmt.Println("debug: deleting atoms not of group")
		// delete all atoms not of the root of either group
		deleteAtomsNotOfGroup(moleculeCopy, []int {group1, group2})

		// check whether each constituent of the double fragment is a hydrocarbon
		isGroup1Hydrocarbon := true
		isGroup2Hydrocarbon := true
		for atomID := range moleculeCopy {
			group := root(moleculeCopy,atomID)
			if moleculeCopy[atomID].element != "H" && moleculeCopy[atomID].element != "C" {
				if group == group1 {
					isGroup1Hydrocarbon = false
				} else if group == group2 {
					isGroup2Hydrocarbon = false
				}
			}
		}


/*
		fmt.Println("debug: replacing HC half")
		// Replace the hydrocarbon half of half-HC fragments with a butyl group
		if isGroup1Hydrocarbon == true && isGroup2Hydrocarbon == false {
			reviseHalfHydrocarbonDoubleFragment(moleculeCopy, borderBonds[i][0], group1)
		} else if isGroup2Hydrocarbon == true && isGroup1Hydrocarbon == false {
			reviseHalfHydrocarbonDoubleFragment(moleculeCopy, borderBonds[i][1], group2)
		}
*/
		dfCopy := copyMolecule(moleculeCopy)
		union(dfCopy,borderBonds[i][0],borderBonds[i][1])

		// save double fragment to array
		dimersSlice[i] = copyMoleculeWithRenumbering(moleculeCopy)
		doubleFragSlice[i] = copyMoleculeWithRenumbering(dfCopy)
		doubleFragIsHydrocarbon[i] = [2]bool {isGroup1Hydrocarbon, isGroup2Hydrocarbon}
	}

	return doubleFragSlice, dimersSlice, doubleFragIsHydrocarbon
}

func deleteAtomsNotOfGroup(atoms map[int]*atom, groups []int) {
	// delete all atoms not of the root of given groups
	atomIDs := make([]int, len(atoms))
	markForDelete := make([]bool, len(atoms))

	i := 0
	for atomID := range atoms {
		atomIDs[i] = atomID
		group := root(atoms,atomID)
		isInValidGroup := false
		for j := 0; j < len(groups); j++ {
			if group == groups[j] {
				isInValidGroup = true
				break
			}
		}
		if !isInValidGroup {
			markForDelete[i] = true
		}
		i++
	}
	for k := 0; k < len(atomIDs); k++ {
		if markForDelete[k] == true {
			delete(atoms, atomIDs[k])
		}
	}
}




// removes bond between two atoms and adds methyls to molecule ends in its place
func removeBondAndCapEnds(atoms map[int]*atom, atom1 int, atom2 int) {
	disconnect(atoms,atom1,atom2)
	capWithMethylGroup(atoms,atom1)
	capWithMethylGroup(atoms,atom2)
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

// currently places atoms w/o coordinates
func capWithMethylGroup(atoms map[int]*atom, atom1 int) {
	// calculate location of new methyl group
	// carbPos := getMethylCoordinates(atoms,atom1,atom2)

	// create new c atom and set parameters
	var carbon atom
	newCarbIndex := len(atoms)+1
	carbon.element = "C"
	carbon.pos = []float64{0.0,0.0,0.0}
	carbon.parent = newCarbIndex
	carbon.atomType = 1
	carbon.bondedAtoms = []int{atom1}
	carbon.treeSize = 1
	carbon.isInFuncGroup = true

	// adjust tree size of root
	atoms[root(atoms,atom1)].treeSize++

	// add to map
	atoms[newCarbIndex] = &carbon
	// bind new carbon to atom1
	atoms[atom1].bondedAtoms = append(atoms[atom1].bondedAtoms, newCarbIndex)
	// join groups
	union(atoms,atom1,newCarbIndex)
	
	for i := 0; i < 3; i++ {
		var hydrogen atom
		newHydrogenIndex := len(atoms)+1
		hydrogen.element = "H"
		hydrogen.pos = []float64{0.0,0.0,0.0}
		hydrogen.parent = newHydrogenIndex
		hydrogen.atomType = 5
		hydrogen.bondedAtoms = []int{newCarbIndex}
		hydrogen.treeSize = 1
		hydrogen.isInFuncGroup = true

		// adjust tree size of root
		atoms[root(atoms,atom1)].treeSize++

		// add to map
		atoms[newHydrogenIndex] = &hydrogen
		// bind new hydrogen to new carbon
		atoms[newCarbIndex].bondedAtoms = append(atoms[newCarbIndex].bondedAtoms, newHydrogenIndex)
		// join groups
		union(atoms,newCarbIndex,newHydrogenIndex)
	}

}


////////////////
// Union Find Alg
////////////////

func union(atoms map[int]*atom, atom1 int, atom2 int) {
	root1 := root(atoms, atom1)
	root2 := root(atoms, atom2)
	if root1 != root2 {
		if atoms[root1].treeSize < atoms[root2].treeSize {
			atoms[root1].parent = root2
			atoms[root2].treeSize += atoms[root1].treeSize
		} else {
			atoms[root2].parent = root1
			atoms[root1].treeSize += atoms[atom1].treeSize
		}
	}
}

func connected(atoms map[int]*atom, atom1 int, atom2 int) bool {
	root1 := root(atoms, atom1)
	root2 := root(atoms, atom2)
	if root1 != root2 {
		return false
	} else {
		return true
	}
}

func validate(atoms map[int]*atom, atom1 int) bool {
	if _, ok := atoms[atom1]; ok {
		return true
	}
	return false
}

func root(atoms map[int]*atom, atom1 int) int {
	isAtomInMap := validate(atoms, atom1)

	if isAtomInMap {
		// store array of visited atoms for path compression afterwards
		var visitedAtoms []int

		// check if atom's parent is equal to atom's parent's parent (i.e. we have reached the top of the tree)
		for atoms[atom1].parent != atoms[atoms[atom1].parent].parent {
			// if not set atom's parent to atom's parent's parent
			atoms[atom1].parent = atoms[atoms[atom1].parent].parent
			// save atom's parent to list of visited atoms to path compress afterwards
			visitedAtoms = append(visitedAtoms, atoms[atom1].parent)
		}

		// compress path
		for _, visitedAtom := range visitedAtoms {
			atoms[visitedAtom].parent = atoms[atom1].parent
		}

		return atoms[atom1].parent
	} else {
		err := errors.New("Call to root(): atom ID key " + strconv.Itoa(atom1) + " was not found in molecule map.")
		log.Fatal(err)
		return -1
	}
}

/////////////////////
// Ring Detection Alg
/////////////////////

var Time = 0

func bridge(atoms map[int]*atom) {

	for atomID, atom := range atoms {
		//fmt.Println(atom.element)
		if atom.visited == false {
			dfs(atoms,atomID)
		}
	}
}

func dfs(atoms map[int]*atom, u int) {

	// mark current node as visited
	atoms[u].visited = true

	// initialize discovery time and low value
	atoms[u].discTime = Time
	atoms[u].minDiscTime = Time
	Time++

	// recurse for all bonded atoms
	for _, v := range atoms[u].bondedAtoms {
		// if bonded atom v is not visited, recurse for it and make it a child of u
		if atoms[v].visited == false {
			atoms[v].parentBF = u
			dfs(atoms, v)

			// check if subtree rooted at v has a connection to an ancestor of u
			atoms[u].minDiscTime = min(atoms[u].minDiscTime, atoms[v].minDiscTime)

			// If the lowest vertex reachable from subtree under v is below u in DFS tree, then u-v is a bridge
			if atoms[v].minDiscTime > atoms[u].discTime {
				// bond = bridge
			} else {
				// bond is not bridge
				atoms[v].isCyclic = true
				atoms[u].isCyclic = true
			}
		} else if v != atoms[u].parentBF { // update min discovery time value of u
			atoms[u].minDiscTime = min(atoms[u].minDiscTime, atoms[v].discTime)
		}
	}
}

