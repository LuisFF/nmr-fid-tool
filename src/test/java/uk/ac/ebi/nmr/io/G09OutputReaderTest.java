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

import java.util.HashMap;

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

            double [][] coupling = new double[14][14];
            coupling[7][8]=-10.3655;
            coupling[7][9]=7.11774;
            coupling[7][10]=2.02412;
            coupling[7][11]=-0.774123;
            coupling[7][12]=0.363124;
            coupling[7][13]=-0.468606;
            coupling[8][9]=13.2608;
            coupling[8][10]=-0.710378;
            coupling[8][11]=0.103987;
            coupling[8][12]= 0.0416588;
            coupling[8][13]=-0.394989;
            coupling[9][10]=5.07155;
            coupling[9][11]=10.6587;
            coupling[9][12]=-0.336642;
            coupling[9][13]=-0.567136;
            coupling[10][11]=-11.1144;
            coupling[10][12]=14.3130;
            coupling[10][13]=-0.139880;
            coupling[11][12]=-0.233422;
            coupling[11][13]=-0.0336012;
            coupling[12][13]=-0.170962;

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

            // simple way to check if structure is well defined
            // it does not work for all molecules
            //TODO make a more reliable structure test
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
            }

            // test the coupling between the hydrogens
            for (int i =7; i< atomContainer.getAtomCount();i++){
                Assert.assertNotNull(atomContainer.getAtom(i).getProperty(G09OutputReader.COUPLING_CONSTANTS));
                for (int j = i+1; j < atomContainer.getAtomCount(); j++){
                    HashMap<IAtom, Double> couplings =(HashMap<IAtom, Double>) atomContainer.getAtom(j)
                            .getProperty(G09OutputReader.COUPLING_CONSTANTS);
                    // get the coupling of atom j with the atom i
                    Assert.assertEquals("Problem reading the proton coupling constantes ",
                            coupling[i][j],couplings.get(atomContainer.getAtom(i)));
//                    System.out.println("Problem reading the proton coupling constantes "+
//                            coupling[i][j] +" " + couplings.get(atomContainer.getAtom(i)).getClass());
                }

            }



            //System.out.println(new SmilesGenerator().createSMILES(atomContainer).toString());
            String smile = "[H]OC(=O)C([H])(N([H])[H])C([H])([H])O[H]";
            Assert.assertEquals("AtomContainer has a different smile",
                    "[H]OC(=O)C([H])(N([H])[H])C([H])([H])O[H]",
                    new SmilesGenerator().createSMILES(atomContainer).toString());

        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }



}
