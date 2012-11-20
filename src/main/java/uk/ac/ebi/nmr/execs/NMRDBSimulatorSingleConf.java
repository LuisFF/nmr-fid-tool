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

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 05/10/2012
 * Time: 16:51
 * To change this template use File | Settings | File Templates.
 */
public class NMRDBSimulatorSingleConf {

    public static void main(String[] args) throws Exception {


        if (args.length < 1)
            throw new Exception("Not enough arguments");

        String filepath = args[0];
        String filename = (new File(filepath)).getName();
        filename = filename.subSequence(0,filename.lastIndexOf(".")).toString();

        // one has to work with the full path otherwise it is not possible to extract the working directory
        String workingDir = (new File((new File(filepath)).getAbsolutePath())).getParent();

        System.out.println("Reading file: " + args[0]);
        System.out.println("working directory" + (new File(workingDir)).getAbsolutePath());
        BufferedReader bufferR = new BufferedReader(new FileReader(new File(filepath)));
        G09OutputReader g09OutputReader = new G09OutputReader(bufferR);
        IChemFile chemFile = new ChemFile();
        chemFile = g09OutputReader.read(chemFile);
        // the reader will read a sequence of atomcontainer and only the last one will have all the information
        IAtomContainer atomContainer = ChemFileManipulator.getAllAtomContainers(chemFile).get(
                ChemFileManipulator.getAllAtomContainers(chemFile).size() - 1);
        IChemModel model = ChemFileManipulator.getAllChemModels(chemFile)
                .get(ChemFileManipulator.getAllChemModels(chemFile).size()-1);
        System.out.println("Energy of the geometry: "+model.getProperty(G09OutputReader.STRUCTURE_ENERGY));
        atomContainer.setProperty(G09OutputReader.STRUCTURE_ENERGY,
                Double.parseDouble((String) model.getProperty(G09OutputReader.STRUCTURE_ENERGY)));
        // calculate the spin systems and corresponding parameters
        SpinSystem spinSystemH = new MagneticPropTool().calculateHSpinSystem(atomContainer);
        // write down the NMR data
        NMRDBSimulatorWriter nmrdbSimulatorWriter = new NMRDBSimulatorWriter(new FileWriter(new File(
                workingDir+"/"+filename+"-simulator.txt"
        )));
        nmrdbSimulatorWriter.write(spinSystemH);



    }
}
