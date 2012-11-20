package uk.ac.ebi.nmr.io;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 27/09/2012
 * Time: 15:30
 * To change this template use File | Settings | File Templates.
 */
public class LasCMDReaderTest {

    /**
     * test if the reader can handle the lascmd input file for a simple molecule like serine
     */
    @Test
    public void testReadSerine() {
        LasCMDReader lasCMDReader = new LasCMDReader(this.getClass()
                .getResourceAsStream("/examples/file_formats/3447.inp"));
        IChemFile chemFile = new ChemFile();


        try {
            chemFile = (IChemFile) lasCMDReader.read(chemFile);
            // the reader will read a sequence of atomcontainer and only the last one will have all the information
            IAtomContainer atomContainer = ChemFileManipulator.getAllAtomContainers(chemFile).get(
                    ChemFileManipulator.getAllAtomContainers(chemFile).size() - 1);

            Assert.assertEquals("Structure was interpreted differently",
                    "[H]C=2N=C1N=C([H])N([H])C1=C(N=2)N([H])[H]",
                    new SmilesGenerator().createSMILES(atomContainer).toString());
            Assert.assertEquals("Structure was interpreted differently",
                    "InChI=1S/C5H5N5/c6-4-3-5(9-1-7-3)10-2-8-4/h1-2H,(H3,6,7,8,9,10)",
                    InChIGeneratorFactory.getInstance().getInChIGenerator(atomContainer).getInchi());

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }
}
