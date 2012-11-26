package uk.ac.ebi.nmr.execs;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.silent.ChemFile;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import uk.ac.ebi.nmr.MagneticPropTool;
import uk.ac.ebi.nmr.SpinSystem;
import uk.ac.ebi.nmr.io.G09OutputReader;
import uk.ac.ebi.nmr.io.L777OutputReader;
import uk.ac.ebi.nmr.io.NMRDBSimulatorWriter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.HashMap;
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
                Pattern REGEX_GAUSSIAN_FILE = Pattern.compile("conf\\_\\d+(\\_solvent)?(-vibcorr)?\\.log");
                return REGEX_GAUSSIAN_FILE.matcher(filename).find();
            }
        };

        File[] gaussianFiles = (new File(workingDir)).listFiles(gaussianFilter);
        if(gaussianFiles.length == 0)
            throw new Exception("No gaussian file found");

        List<IAtomContainer> vacuumSimulations = new ArrayList<IAtomContainer>();
        List<IAtomContainer> solventSimulations = new ArrayList<IAtomContainer>();
        List<IAtomContainer> vacuumVibcorrSimulations = new ArrayList<IAtomContainer>();
        List<IAtomContainer> solventVibcorrSimulations = new ArrayList<IAtomContainer>();
        System.out.println(gaussianFiles.length);
        for (File file : gaussianFiles) {
            BufferedReader bufferR = new BufferedReader(new FileReader(file));
            G09OutputReader g09OutputReader = new G09OutputReader(bufferR);
            L777OutputReader l777OutputReader = new L777OutputReader(bufferR);
            IChemFile chemFile = new ChemFile();
            if(file.getName().contains("vibcorr")){
                chemFile = l777OutputReader.read(chemFile);
            }else{
                chemFile = g09OutputReader.read(chemFile);
            }
            // the reader will read a sequence of atomcontainer and only the last one will have all the information
            IChemModel model = ChemFileManipulator.getAllChemModels(chemFile)
                    .get(ChemFileManipulator.getAllChemModels(chemFile).size() - 1);

            if(file.getName().contains("vibcorr")){
                model.setProperty(G09OutputReader.NORMAL_TERMINATION,
                        (Boolean) model.getProperty(L777OutputReader.NORMAL_TERMINATION));
            }

            System.out.println("processing file: "+file.getName());
            System.out.println("Normal termination: "+
                    ((Boolean) model.getProperty(G09OutputReader.NORMAL_TERMINATION)));

            if ((Boolean) model.getProperty(G09OutputReader.NORMAL_TERMINATION)) {
                IAtomContainer atomContainer = ChemFileManipulator.getAllAtomContainers(chemFile)
                        .get(ChemFileManipulator.getAllAtomContainers(chemFile).size() - 1);
                if(file.getName().contains("vibcorr")){
                    for(IAtom atom : atomContainer.atoms()){
                        atom.setProperty(G09OutputReader.MAGNETIC_TENSOR,
                                (Double) atom.getProperty(L777OutputReader.MAGNETIC_TENSOR));
                        atom.setProperty(G09OutputReader.COUPLING_CONSTANTS,
                                (HashMap<IAtom,Double>) atom.getProperty(L777OutputReader.COUPLING_CONSTANTS));
                    }
                    System.out.println(model.getProperty(L777OutputReader.STRUCTURE_ENERGY));
                    atomContainer.setProperty(G09OutputReader.STRUCTURE_ENERGY,
                            Double.parseDouble((String) model.getProperty(L777OutputReader.STRUCTURE_ENERGY)));
                } else {
                    atomContainer.setProperty(G09OutputReader.STRUCTURE_ENERGY,
                            Double.parseDouble((String) model.getProperty(G09OutputReader.STRUCTURE_ENERGY)));
                }
                System.out.println("Energy of the geometry: " +
                        atomContainer.getProperty(G09OutputReader.STRUCTURE_ENERGY));
                //sort the atomcontainer
                if (file.getName().contains("solvent")) {
                    if(file.getName().contains("vibcorr")){
                        solventVibcorrSimulations.add(atomContainer);
                    } else {
                        solventSimulations.add(atomContainer);
                    }
                } else {
                    if(file.getName().contains("vibcorr")){
                        vacuumVibcorrSimulations.add(atomContainer);
                    } else {
                        vacuumSimulations.add(atomContainer);
                    }
                }
            }
            // write down the NMR data
        }
        if(vacuumSimulations.size()==0)
            throw new Exception("No simulations in vacuum conditions found");
        if(solventSimulations.size()==0)
            throw new Exception("No simulations in solvent conditions found");
        if(solventVibcorrSimulations.size()==0)
            throw new Exception("No simulations in solvent and vibrational corrections conditions found");
        if(vacuumVibcorrSimulations.size()==0)
            throw new Exception("No simulations in vacuum and vibrational corrections conditions found");

        SpinSystem averageSpinSystemVacuum = new MagneticPropTool().calculateAverageHSpinSystem(vacuumSimulations);
        SpinSystem averageSpinSystemSolvent = new MagneticPropTool().calculateAverageHSpinSystem(solventSimulations);
        SpinSystem averageSpinSystemVacuumVibbcorr = new MagneticPropTool()
                .calculateAverageHSpinSystem(vacuumVibcorrSimulations);
        SpinSystem averageSpinSystemSolventVibbcorr = new MagneticPropTool()
                .calculateAverageHSpinSystem(solventVibcorrSimulations);

        NMRDBSimulatorWriter nmrdbSimulatorWriter = new NMRDBSimulatorWriter(new FileWriter(new File(
                workingDir+"/avgSpectra-vacuum-sim.txt")));
        nmrdbSimulatorWriter.write(averageSpinSystemVacuum);

        nmrdbSimulatorWriter = new NMRDBSimulatorWriter(new FileWriter(new File(
                workingDir+"/avgSpectra-vacuum-vibcorr-sim.txt")));
        nmrdbSimulatorWriter.write(averageSpinSystemVacuumVibbcorr);

        nmrdbSimulatorWriter = new NMRDBSimulatorWriter(new FileWriter(new File(
                workingDir+"/avgSpectra-solvent-sim.txt")));
        nmrdbSimulatorWriter.write(averageSpinSystemSolvent);

        nmrdbSimulatorWriter = new NMRDBSimulatorWriter(new FileWriter(new File(
                workingDir+"/avgSpectra-solvent-vibcorr-sim.txt")));
        nmrdbSimulatorWriter.write(averageSpinSystemSolventVibbcorr);

    }
}
