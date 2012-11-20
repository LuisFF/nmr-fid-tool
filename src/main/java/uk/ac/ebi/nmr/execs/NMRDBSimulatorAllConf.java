package uk.ac.ebi.nmr.execs;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.silent.ChemFile;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import uk.ac.ebi.nmr.MagneticPropTool;
import uk.ac.ebi.nmr.SpinSystem;
import uk.ac.ebi.nmr.io.G09OutputReader;
import uk.ac.ebi.nmr.io.NMRDBSimulatorWriter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 09/10/2012
 * Time: 14:43
 * To change this template use File | Settings | File Templates.
 */
public class NMRDBSimulatorAllConf {


    public static void main(String[] args) throws Exception {


        if (args.length < 1)
            throw new Exception("Not enough arguments");
        if (!(new File(args[0]).isDirectory()))
            throw new Exception("The argument is not a directory");

        String workingDir = (new File(args[0])).getAbsolutePath();


        System.out.println("working directory: " + workingDir);

        FilenameFilter gaussianFilter = new FilenameFilter() {
            @Override
            public boolean accept(File file, String filename) {
                Pattern REGEX_GAUSSIAN_FILE = Pattern.compile("conf\\_\\d+(\\_solvent)?\\.log");
                return REGEX_GAUSSIAN_FILE.matcher(filename).find();
            }
        };

        File[] gaussianFiles = (new File(workingDir)).listFiles(gaussianFilter);
        if(gaussianFiles.length == 0)
            throw new Exception("No gaussian file found");

        List<IAtomContainer> vacuumSimulations = new ArrayList<IAtomContainer>();
        List<IAtomContainer> solventSimulations = new ArrayList<IAtomContainer>();
        System.out.println(gaussianFiles.length);
        for (File file : gaussianFiles) {
            BufferedReader bufferR = new BufferedReader(new FileReader(file));
            G09OutputReader g09OutputReader = new G09OutputReader(bufferR);
            IChemFile chemFile = new ChemFile();
            chemFile = g09OutputReader.read(chemFile);
            // the reader will read a sequence of atomcontainer and only the last one will have all the information
            IChemModel model = ChemFileManipulator.getAllChemModels(chemFile)
                    .get(ChemFileManipulator.getAllChemModels(chemFile).size() - 1);
            System.out.println((Boolean) model.getProperty(G09OutputReader.NORMAL_TERMINATION));
            if ((Boolean) model.getProperty(G09OutputReader.NORMAL_TERMINATION)) {
                System.out.println("processing file: "+file.getName());
                IAtomContainer atomContainer = ChemFileManipulator.getAllAtomContainers(chemFile)
                        .get(ChemFileManipulator.getAllAtomContainers(chemFile).size() - 1);
                System.out.println("Energy of the geometry: " + model.getProperty(G09OutputReader.STRUCTURE_ENERGY));
                atomContainer.setProperty(G09OutputReader.STRUCTURE_ENERGY,
                        Double.parseDouble((String) model.getProperty(G09OutputReader.STRUCTURE_ENERGY)));
                //sort the atomcontainer
                if (file.getName().contains("solvent")) {
                    solventSimulations.add(atomContainer);
                } else {
                    vacuumSimulations.add(atomContainer);
                }
            }
            // write down the NMR data
        }
        if(vacuumSimulations.size()==0)
            throw new Exception("No simulations in vacuum conditions found");
        if(solventSimulations.size()==0)
            throw new Exception("No simulations in solvent conditions found");

        SpinSystem averageSpinSystemVacuum = new MagneticPropTool().calculateAverageHSpinSystem(vacuumSimulations);
        SpinSystem averageSpinSystemSolvent = new MagneticPropTool().calculateAverageHSpinSystem(solventSimulations);

        NMRDBSimulatorWriter nmrdbSimulatorWriter = new NMRDBSimulatorWriter(new FileWriter(new File(
                workingDir+"/averageSpectra-vacuum-simulator.txt")));
        nmrdbSimulatorWriter.write(averageSpinSystemVacuum);

        nmrdbSimulatorWriter = new NMRDBSimulatorWriter(new FileWriter(new File(
                workingDir+"/averageSpectra-solvent-simulator.txt")));
        nmrdbSimulatorWriter.write(averageSpinSystemSolvent);
    }
}
