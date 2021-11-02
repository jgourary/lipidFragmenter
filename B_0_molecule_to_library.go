package main

import (
	"fmt"
	"io/ioutil"
	"log"
	"os"
	filepath2 "path/filepath"
	"strconv"
)


func generateAtomCodeDictFile(dir string, outFile string, outName string) {
	fileInfo, err := ioutil.ReadDir(dir)
	if err != nil {
		fmt.Println("failed to read directory: " + dir)
		log.Fatal(err)
	}
	_ = os.MkdirAll(dir, 0755)
	thisFile, err := os.Create(filepath2.Join(outFile, outName))
	if err != nil {
		fmt.Println("Failed to create new atom code file: " + outFile)
		log.Fatal(err)
	}


	for i := 0; i < len(fileInfo); i++ {
		if filepath2.Ext(fileInfo[i].Name()) == ".txyz" {
			txyzFilePath := filepath2.Join(dir, fileInfo[i].Name())
			molName, atoms := loadLipid(txyzFilePath)
			atomCodeToTypeMap := getAtomCodeToAtomTypeMap(atoms)
			for k, v := range atomCodeToTypeMap {
				_, _ = thisFile.WriteString(k + "\t" + strconv.Itoa(v) + "\t" + molName + "\n")
			}
		}
	}
}


func getAtomCodeToAtomTypeMap(atoms map[int]*atom) map[string]int {
	atomCodeToTypeMap := make(map[string]int)
	for atomID := range atoms {
		atomCodeToTypeMap[getAtomCode(atoms, atomID)] = atoms[atomID].atomType
	}
	return atomCodeToTypeMap
}

