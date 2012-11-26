package uk.ac.ebi.nmr.execs;

import org.openscience.cdk.exception.CDKException;
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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 05/10/2012
 * Time: 16:51
 * To change this template use File | Settings | File Templates.
 */
public class NMRDBSimulatorSingleConf {

    private final static Pattern REGEX_VIBRATIONAL_CORRECTION_INPUT_FILE = Pattern.compile("vibcorr");

    public static void main(String[] args) throws Exception {


        if (args.length < 1)
            throw new Exception("Not enough arguments");

        String filepath = args[0];
        String filename = (new File(filepath)).getName();
        filename = filename.subSequence(0,filename.lastIndexOf(".")).toString();

        // one has to work with the full path otherwise it is not possible to extract the working directory
        String workingDir = (new File((new File(filepath)).getAbsolutePath())).getParent();

        System.out.println("Reading file: " + args[0]);
        System.out.println("working directory: " + (new File(workingDir)).getAbsolutePath());
        BufferedReader bufferR = new BufferedReader(new FileReader(new File(filepath)));

        G09OutputReader g09OutputReader = new G09OutputReader(bufferR);
        L777OutputReader l777OutputReader = new L777OutputReader(bufferR);
        IChemFile chemFile = new ChemFile();
        if(!REGEX_VIBRATIONAL_CORRECTION_INPUT_FILE.matcher(filename).find()){
            chemFile = g09OutputReader.read(chemFile);
        }else{
            chemFile = l777OutputReader.read(chemFile);
        }
        // the reader will read a sequence of atomcontainer and only the last one will have all the information
        IAtomContainer atomContainer = ChemFileManipulator.getAllAtomContainers(chemFile).get(
                ChemFileManipulator.getAllAtomContainers(chemFile).size() - 1);
        IChemModel model = ChemFileManipulator.getAllChemModels(chemFile)
                .get(ChemFileManipulator.getAllChemModels(chemFile).size()-1);

        // set G09OutputReader properties if reader is L777OutputReader
        // these properties are required by the MagneticPropTool
        if(REGEX_VIBRATIONAL_CORRECTION_INPUT_FILE.matcher(filename).find()){
            for(IAtom atom : atomContainer.atoms()){
                atom.setProperty(G09OutputReader.MAGNETIC_TENSOR,
                    (Double) atom.getProperty(L777OutputReader.MAGNETIC_TENSOR));
                atom.setProperty(G09OutputReader.COUPLING_CONSTANTS,
                        (HashMap<IAtom,Double>) atom.getProperty(L777OutputReader.COUPLING_CONSTANTS));
            }
            atomContainer.setProperty(G09OutputReader.STRUCTURE_ENERGY,
                    Double.parseDouble((String) model.getProperty(L777OutputReader.STRUCTURE_ENERGY)));
        } else {
            atomContainer.setProperty(G09OutputReader.STRUCTURE_ENERGY,
                    Double.parseDouble((String) model.getProperty(G09OutputReader.STRUCTURE_ENERGY)));
        }

        System.out.println("Energy of the geometry: "+atomContainer.getProperty(G09OutputReader.STRUCTURE_ENERGY));


        //write the gaussian file for the vibrational corrections
        if(!REGEX_VIBRATIONAL_CORRECTION_INPUT_FILE.matcher(filename).find())
            writeGaussian(filename, workingDir, atomContainer);
        // calculate the spin systems and corresponding parameters
        SpinSystem spinSystemH = new MagneticPropTool().calculateHSpinSystem(atomContainer);
        // write down the NMR data
        NMRDBSimulatorWriter nmrdbSimulatorWriter = new NMRDBSimulatorWriter(new FileWriter(new File(
                workingDir+"/"+filename+"-simulator.txt"
        )));
        nmrdbSimulatorWriter.write(spinSystemH);

    }
    private static void writeGaussian(String filename, String workingDir,
                                      IAtomContainer atomContainer) throws IOException, CDKException {
            //
            String header =
                    "! Iop(7/95=298150)  ====> T=298150/1000\n" +
                    "! IOp(7/96=30)      ====> Delta = 30/1000\n" +
                    "! IOp(7/97=1)       ====> Do spin-spin (0 for shielding only)\n"+
                    "%mem=40MW\n" +
                    "%nproc=12\n" +
                    "%chk="+ filename+ "-vibcorr.chk\n" +
                    "%subst l777 .\n" +
                    "# B3LYP/6-311++G** Freq(Anha) \n" +
                    "# extralink=L777 Iop(7/95=298150) IOp(7/96=30) IOp(7/97=1)\n\n" +
                    "Geometry optimized by Gaussian = " +
                    atomContainer.getProperty(G09OutputReader.STRUCTURE_ENERGY)+
                    "\n[vibrational corrections]\n" +
                    "\n" +
                    "0,1\n";
            FileWriter fstream = new FileWriter(workingDir + "/" + filename + "-vibcorr.inp");
            BufferedWriter out = new BufferedWriter(fstream);
            out.write(header);
            String atomCoordinates = "";
            for (IAtom atom : atomContainer.atoms()) {
                atomCoordinates += atom.getSymbol() + "," +
                        atom.getPoint3d().x + "," +
                        atom.getPoint3d().y + "," +
                        atom.getPoint3d().z + "\n";
            }
            out.write(atomCoordinates);
            out.write("\n\n");
            out.close();
    }
}
