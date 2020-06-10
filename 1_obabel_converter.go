package main

import (
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	filepath2 "path/filepath"
	"strings"
	"sync"
)

func obabelConversion(directory string, ext1 string, ext2 string, addHydrogens bool) {
	// Read in all files in dir
	fileInfo, err := ioutil.ReadDir(directory)
	if err != nil {
		fmt.Println("failed to read directory: " + directory)
		log.Fatal(err)
	}


	// Iterate through all items in directory
	frame := []int{0,goRoutinesBatchSize-1}

	for frame[0] < len(fileInfo) {

		wg := sync.WaitGroup{}

		for i := frame[0]; i < frame[1]; i++ {
			// if item is a Dir (as it should be unless the end user tampered with the directory manually...)
			if fileInfo[i].IsDir() {
				subFileInfo, err := ioutil.ReadDir(filepath2.Join(directory, fileInfo[i].Name()))
				if err != nil {
					fmt.Println("failed to read directory: " + directory)
					log.Fatal(err)
				}
				for j := 0; j < len(subFileInfo); j++ {
					if filepath2.Ext(subFileInfo[j].Name()) == ext2 {
						os.Remove(filepath2.Join(directory, fileInfo[i].Name(), subFileInfo[j].Name()))
					} else if filepath2.Ext(subFileInfo[j].Name()) == ext1 {
						baseName := strings.Split(subFileInfo[j].Name(), ".")[0]
						convName := baseName + ext2
						basePath := filepath2.Join(directory, fileInfo[i].Name(), subFileInfo[j].Name())
						convPath := filepath2.Join(directory, fileInfo[i].Name(), convName)

						wg.Add(1)
						go obabelWrapper(basePath, convPath, addHydrogens, &wg)

					}
				}
			}
		}
		frame[0] += goRoutinesBatchSize
		frame[1] += goRoutinesBatchSize
		frame[1] = min(frame[1], len(fileInfo))
		wg.Wait()
	}

}

func obabelWrapper(path1 string, path2 string, addHydrogens bool, wg *sync.WaitGroup) {
	obabel := "C:\\Program Files (x86)\\OpenBabel-2.3.1\\obabel.exe"
	var cmdArgs []string
	if addHydrogens == true {
		cmdArgs = []string{path1, "-O", path2, "-h"}
	} else {
		cmdArgs = []string{path1, "-O", path2}
	}
	//fmt.Println(path2)
	//cmdstring := obabel + " -i " + path1 + " -o " + path2
	out, err := exec.Command(obabel, cmdArgs...).CombinedOutput()
	//fmt.Println(string(out))
	wg.Done()
	if err != nil {
		fmt.Println(string(out))
		fmt.Println(err)
		log.Fatal(err)
	}
}
