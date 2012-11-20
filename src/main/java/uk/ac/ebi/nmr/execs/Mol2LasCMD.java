package uk.ac.ebi.nmr.execs;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import uk.ac.ebi.nmr.io.LasCMDWriter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 01/08/2012
 * Time: 16:45
 * To change this template use File | Settings | File Templates.
 */
public class Mol2LasCMD {
    public static void main(String[] args) throws Exception {
        if (args.length < 1)
            throw new Exception("Not enough arguments");

        String filepath = args[0];
        String workingDir = (new File(filepath)).getParent();
        String filename = (new File(filepath)).getName().split("\\.")[0];



        BufferedReader bufferR = new BufferedReader(new FileReader(new File(filepath)));
        MDLV2000Reader reader = new MDLV2000Reader(bufferR);
        IChemFile chemFile = new ChemFile();
        chemFile = (IChemFile) reader.read(chemFile);

        LasCMDWriter writer = new LasCMDWriter(new FileWriter(new File(workingDir+"/"+filename+".inp")));
        writer.write(ChemFileManipulator.getAllAtomContainers(chemFile).get(0));
    }
}
