package uk.ac.ebi.nmr.io;

import junit.framework.Assert;
import org.junit.Test;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.alignment.KabschAlignment;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 13/08/2012
 * Time: 11:43
 * To change this template use File | Settings | File Templates.
 */
public class G09OutputReaderTest {

    @Test
    public void testReadSerine(){
        try {
            double [] chemShifts = {229.1109, 122.8042, 1.9914, -85.3661,
                    115.4538, 292.7408, 121.7134, 31.0980, 30.6853,
                    28.5072, 28.1297, 28.7746, 29.2412, 26.1421};

            //get the data from the gaussian file
//            BufferedReader bufferR = new BufferedReader(new FileReader(
//                    new File("/Users/ldpf/NetBeansProjects/cdk-playground/src/test/java/resources/examples/file_formats/conf_0.log")));
////                    new File("/Volumes/ldpf/quantum_chemistry/compounds/adenine.log")));
////                    new File("/Users/ldpf/Data/NMR/QM/Merche-Simulations/adenine.log")));
////                    new File("/Users/ldpf/Downloads/adenine.log")));
////                    new File("/Volumes/ldpf/quantum_chemistry/compounds/14hydrthoic.log")));
//            G09OutputReader g09OutputReader = new G09OutputReader(bufferR);

            G09OutputReader g09OutputReader = new G09OutputReader(this.getClass()
                    .getResourceAsStream("/examples/file_formats/conf_0.log"));
            IChemFile chemFile = new ChemFile();
            chemFile = (IChemFile) g09OutputReader.read(chemFile);
            // the reader will read a sequence of atomcontainer and only the last one will have all the information
            IAtomContainer atomContainer = ChemFileManipulator.getAllAtomContainers(chemFile).get(
                        ChemFileManipulator.getAllAtomContainers(chemFile).size()-1);

            Assert.assertEquals("Structure was interpreted differently",
                    "[H]OC(=O)C([H])(N([H])[H])C([H])([H])O[H]",
                    new SmilesGenerator().createSMILES(atomContainer).toString());

            // debug, check if geometry changed
                IAtomContainer container = ChemFileManipulator.getAllAtomContainers(chemFile).get(0);
                IAtomContainer container1 = ChemFileManipulator.getAllAtomContainers(chemFile)
                        .get(ChemFileManipulator.getAllAtomContainers(chemFile).size()-1);
                KabschAlignment kabschAlignment = new KabschAlignment(container,container1);
                kabschAlignment.align();
            Assert.assertTrue("There is no diference between first and last geometry",kabschAlignment.getRMSD()>0);
            // debug check if the values were correctly read
            for (int i = 0 ; i < atomContainer.getAtomCount();i++){
                IAtom atom = atomContainer.getAtom(i);
                Assert.assertEquals("Problem reading the magnetic tensor values", chemShifts[i],
                        (Double) atom.getProperty(G09OutputReader.MAGNETIC_TENSOR));
                if(atom.getProperty(G09OutputReader.COUPLING_CONSTANTS) != null)
                System.out.println(atom.getProperty(G09OutputReader.COUPLING_CONSTANTS).toString());
            }



            System.out.println(new SmilesGenerator().createSMILES(atomContainer).toString());
            String smile = "[H]OC(=O)C([H])(N([H])[H])C([H])([H])O[H]";
            Assert.assertEquals("AtomContainer has a different smile",
                    "[H]OC(=O)C([H])(N([H])[H])C([H])([H])O[H]",
                    new SmilesGenerator().createSMILES(atomContainer).toString());

        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }


}
