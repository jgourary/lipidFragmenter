package main

import (
	"errors"
	"log"
	"math"
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
}

///////////////////
// Functions
///////////////////

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
// groupSlice stores the original atomIDs of all atoms in the group
func groupMap(atoms map[int]*atom, group []int) (map[int]int, []int) {
	groupMap := make(map[int]int)
	groupSlice := make([]int,len(group))
	for i, atomID := range group {
		groupMap[atomID] = i+1
		groupSlice[i] = atomID
	}

	return groupMap, groupSlice
}

// removes bond between two atoms and adds hydrogens to molecule ends in its place
func removeBond(atoms map[int]*atom, atom1 int, atom2 int) {
	disconnect(atoms,atom1,atom2)
	capWithMethylGroup(atoms,atom1,atom2)
	capWithMethylGroup(atoms,atom2,atom1)
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

// currently only adds carbon atom - does not add accompanying hydrogens
func capWithMethylGroup(atoms map[int]*atom, atom1 int, atom2 int) {
	// calculate location of new h atom
	carbPos := getCarbonCoordinates(atoms,atom1,atom2)

	// create new h atom and set parameters
	var carbon atom
	index := len(atoms)+1
	carbon.element = "C"
	carbon.pos = carbPos
	carbon.parent = index
	carbon.atomType = 1
	carbon.bondedAtoms = []int{atom1}
	carbon.treeSize = 1
	carbon.isInFuncGroup = true

	// adjust tree size of root
	atoms[root(atoms,atom1)].treeSize++

	// add to map

	atoms[index] = &carbon
	// bind new carbon to atom1
	atoms[atom1].bondedAtoms = append(atoms[atom1].bondedAtoms, index)
	// join groups
	union(atoms,atom1,index)
}

func getCarbonCoordinates(atoms map[int]*atom, atom1 int, atom2 int) []float64 {
	atom1coords := atoms[atom1].pos
	atom2coords := atoms[atom2].pos

	distance := 0.0
	for i := 0; i < 3; i++ {
		distance += math.Pow(atom2coords[i] - atom1coords[i],2)
	}
	distance = math.Sqrt(distance)

	r := make([]float64,3)
	for i:=0; i < 3; i++ {
		r[i] = (carbonCarbonBondDistance * (atom2coords[i] - atom1coords[i]) / distance) + atom1coords[i]
	}

	return r
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
	} else {
		newErr := errors.New("Atom key " + strconv.Itoa(atom1) + " was not found in map")
		log.Fatal(newErr)
	}
	return false
}

func root(atoms map[int]*atom, atom1 int) int {
	validate(atoms, atom1)
	// store array of visited atoms for path compression afterwards
	var visitedAtoms []int

	// check if atom's parent is equal to atom's parent's parent (i.e. we have reached the top of the tree)
	for atoms[atom1].parent != atoms[atoms[atom1].parent].parent {
		// if not set atom's parent to atom's parent's parent
		atoms[atom1].parent = atoms[atoms[atom1].parent].parent
		// save atom's parent to list of visited atoms to path compress afterwards
		visitedAtoms = append(visitedAtoms,atoms[atom1].parent)
	}

	// compress path
	for _,visitedAtom := range visitedAtoms {
		atoms[visitedAtom].parent = atoms[atom1].parent
	}

	return atoms[atom1].parent
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

