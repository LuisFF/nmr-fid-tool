package uk.ac.ebi.nmr.io;

import org.junit.Test;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.io.MDLV2000Reader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;


/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 01/08/2012
 * Time: 10:07
 * To change this template use File | Settings | File Templates.
 */
public class LasCMDWriterTest {

    @Test
    public void writerTest(){
        try {
            BufferedReader bufferR = new BufferedReader(new FileReader(
                    new File("/Users/ldpf/NetBeansProjects/cdk-playground/src/test/java/resources/3343.mol")));
            MDLV2000Reader reader = new MDLV2000Reader(bufferR);
            IChemFile chemFile = new ChemFile();
            chemFile = (IChemFile) reader.read(chemFile);
//            LasCMDWriter writer = new LasCMDWriter(new FileWriter(new File("/Users/ldpf/NetBeansProjects/cdk-playground/src/test/java/resources/3343_test.inp")));
//            writer.write(ChemFileManipulator.getAllAtomContainers(chemFile).get(0));
//            System.out.println(chemFile.getID());
//            ChemFileManipulator.getAllAtomContainers(chemFile);
//            LasCMDWriter writer = new LasCMDWriter(container);
//            writer.write();


        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }
}
