import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.silent.ChemFile;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import uk.ac.ebi.nmr.MagneticPropTool;
import uk.ac.ebi.nmr.SpinSystem;
import uk.ac.ebi.nmr.io.G09OutputReader;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 28/08/2012
 * Time: 11:35
 * To change this template use File | Settings | File Templates.
 */
public class MagneticPropToolTest {

    IAtomContainer atomContainer;

    @Before
    public void initialize() {

        BufferedReader bufferR = new BufferedReader(
                new InputStreamReader(this.getClass()
                        .getResourceAsStream("/examples/file_formats/conf_0.log")));
//                    new FileReader(
//                new File("/Users/ldpf/NetBeansProjects/cdk-playground/src/test/java/resources/examples/file_formats/conf_0.log")));
//                    new File("/Volumes/ldpf/quantum_chemistry/compounds/adenine.log")));
//                    new File("/Users/ldpf/Data/NMR/QM/Merche-Simulations/adenine.log")));
//                    new File("/Users/ldpf/Downloads/adenine.log")));
//                    new File("/Volumes/ldpf/quantum_chemistry/compounds/14hydrthoic.log")));
//                            new File("/Volumes/ldpf/quantum_chemistry/compounds/acetovanillone.log")));


        try {

            G09OutputReader g09OutputReader = new G09OutputReader(bufferR);
            IChemFile chemFile = new ChemFile();
            chemFile = g09OutputReader.read(chemFile);
            IChemModel model = ChemFileManipulator.getAllChemModels(chemFile).get(
                    ChemFileManipulator.getAllChemModels(chemFile).size() - 1);
            // the reader will read a sequence of atomcontainer and only the last one will have all the information
            atomContainer = ChemFileManipulator.getAllAtomContainers(chemFile).get(
                    ChemFileManipulator.getAllAtomContainers(chemFile).size() - 1);
            // at the moment I have to set the energy property manually in the atomContainer obtained from the model
            atomContainer.setProperty(G09OutputReader.STRUCTURE_ENERGY,
                    Double.parseDouble((String) model.getProperty(G09OutputReader.STRUCTURE_ENERGY)));
            g09OutputReader.close();

        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e1) {
            e1.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }


    @Test
    public void testEquivClassPartitioning() {
        try {
            SpinSystem spinSystem = new MagneticPropTool().calculateHSpinSystem(atomContainer);
            int count = 0;
            HashSet<Integer> classes = new HashSet<Integer>();

            for (IAtom atom : spinSystem.getAtomContainer().atoms()) {
                if (atom.getSymbol().equals("C")) {
                    Assert.assertNotNull("No partitioning done", atom.getProperty(MagneticPropTool.EQUIVALENT_CLASS));
                    if (!classes.contains((Integer) atom.getProperty(MagneticPropTool.EQUIVALENT_CLASS))) {
                        count++;
                        classes.add((Integer) atom.getProperty(MagneticPropTool.EQUIVALENT_CLASS));
                    }
                }
            }
            Assert.assertEquals("Problem with partitioning", 3, count);
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }

    @Test
    public void testHSpinSystem() {
        try {

            SpinSystem spinSystemH = new MagneticPropTool().calculateHSpinSystem(atomContainer);
            Assert.assertNotNull("Problem creating the spin system", spinSystemH);

            HashMap<Integer, Double> magneticTensors = new HashMap<Integer, Double>();
            magneticTensors.put(1, 28.5072);
            magneticTensors.put(2, 28.45215);
            magneticTensors.put(0, 0.0);

            double couplingConstant = 7.9;

            // check if the partitioning went correct
            Assert.assertEquals("Different grouping of spin systems", 3, spinSystemH.getNumberOfSpins().length);
            boolean numberOfSpinsBiggerThanZero = false;
            // check if one has the same distribution of spins per spin group
            for (int i : spinSystemH.getNumberOfSpins()) {
                Assert.assertTrue("Number of spins per spin group is wrong", magneticTensors.containsKey(i));
                numberOfSpinsBiggerThanZero |= (i>0);
            }
            Assert.assertTrue("Problem creating the sim system", numberOfSpinsBiggerThanZero);
            // check if the magnetic tensor is assigned to the correct group
            for (int i = 0; i < spinSystemH.getNbOfChemEquivalentSpinGroups(); i++) {
                Assert.assertEquals("Wrong magnetic tensor assignment ",
                        magneticTensors.get(spinSystemH.getNumberOfSpins()[i]),
                        spinSystemH.getMagneticTensors()[i]);
                // check if the coupling constant is correctly calculated
                for (int j = 0; j < spinSystemH.getNbOfChemEquivalentSpinGroups(); j++) {
                    if (i != j && spinSystemH.getNumberOfSpins()[i] != 0 && spinSystemH.getNumberOfSpins()[j] != 0) {

                        if (spinSystemH.getConstantJMatrix()[i][j] != 0)
                            Assert.assertEquals("Wrong coupling constant", 7.865124940872192,
                                    spinSystemH.getConstantJMatrix()[i][j]);
                    }

                }
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }


    @Test
    public void testAverageHSpinSystem() {

        try {

            IAtomContainer atomContainer1 = SilentChemObjectBuilder.getInstance()
                    .newInstance(IAtomContainer.class, atomContainer.clone());
            IAtomContainer atomContainer2 = SilentChemObjectBuilder.getInstance()
                    .newInstance(IAtomContainer.class, atomContainer.clone());
            // set all properties that were not copied by the cloning of the atomContainer

            atomContainer1.setProperty(G09OutputReader.STRUCTURE_ENERGY,
                    (Double) atomContainer.getProperty(G09OutputReader.STRUCTURE_ENERGY)+0.000001);
            atomContainer2.setProperty(G09OutputReader.STRUCTURE_ENERGY,
                    (Double) atomContainer.getProperty(G09OutputReader.STRUCTURE_ENERGY)+0.05);
            for (int i = 0 ; i < atomContainer.getAtomCount(); i++){
                if(atomContainer.getAtom(i).getProperty(G09OutputReader.MAGNETIC_TENSOR)!= null) {
                    atomContainer1.getAtom(i).setProperty(G09OutputReader.MAGNETIC_TENSOR,
                            atomContainer.getAtom(i).getProperty(G09OutputReader.MAGNETIC_TENSOR));
                    atomContainer2.getAtom(i).setProperty(G09OutputReader.MAGNETIC_TENSOR,
                            atomContainer.getAtom(i).getProperty(G09OutputReader.MAGNETIC_TENSOR));
                }
                if(atomContainer.getAtom(i).getProperty(G09OutputReader.COUPLING_CONSTANTS)!=null){
                    HashMap<IAtom, Float> coupledAtoms = (HashMap<IAtom, Float>) atomContainer
                            .getAtom(i).getProperty(G09OutputReader.COUPLING_CONSTANTS);
                    for(IAtom atom : coupledAtoms.keySet()){
                        if(atomContainer1.getAtom(i).getProperty(G09OutputReader.COUPLING_CONSTANTS) !=null){
                            atomContainer1.getAtom(i)
                                    .setProperty(G09OutputReader.COUPLING_CONSTANTS, new HashMap<IAtom, Float>());
                            atomContainer2.getAtom(i)
                                    .setProperty(G09OutputReader.COUPLING_CONSTANTS, new HashMap<IAtom, Float>());
                        }
                        ((HashMap<IAtom, Float>) atomContainer1.getAtom(i)
                                .getProperty(G09OutputReader.COUPLING_CONSTANTS)).put(atomContainer1.getAtom(
                                atomContainer.getAtomNumber(atom)),
                                coupledAtoms.get(atom));
                        ((HashMap<IAtom, Float>) atomContainer2.getAtom(i)
                                .getProperty(G09OutputReader.COUPLING_CONSTANTS)).put(atomContainer2.getAtom(
                                atomContainer.getAtomNumber(atom)),
                                coupledAtoms.get(atom));
                    }
                }
            }



            for (int i = 0 ; i < atomContainer.getAtomCount(); i++){
                if(atomContainer.getAtom(i).getSymbol().equals("H")){
//                    System.out.println("Magnetic tensors before: "+
//                            atomContainer.getAtom(i).getProperty(G09OutputReader.MAGNETIC_TENSOR)+ " "+
//                            atomContainer1.getAtom(i).getProperty(G09OutputReader.MAGNETIC_TENSOR)+ " "+
//                            atomContainer2.getAtom(i).getProperty(G09OutputReader.MAGNETIC_TENSOR));
                    atomContainer1.getAtom(i).setProperty(G09OutputReader.MAGNETIC_TENSOR,
                            (Double) atomContainer1.getAtom(i).getProperty(G09OutputReader.MAGNETIC_TENSOR)+0.5);
                    atomContainer2.getAtom(i).setProperty(G09OutputReader.MAGNETIC_TENSOR,
                            (Double) atomContainer2.getAtom(i).getProperty(G09OutputReader.MAGNETIC_TENSOR)+0.3);
//                    System.out.println("Magnetic tensors after: "+
//                            atomContainer.getAtom(i).getProperty(G09OutputReader.MAGNETIC_TENSOR)+ " "+
//                            atomContainer1.getAtom(i).getProperty(G09OutputReader.MAGNETIC_TENSOR)+ " "+
//                            atomContainer2.getAtom(i).getProperty(G09OutputReader.MAGNETIC_TENSOR));
                }
                // iterate through the coupled atoms...
            }
            List<IAtomContainer> containerList = new ArrayList<IAtomContainer>();

            containerList.add(atomContainer);
            containerList.add(atomContainer1);
            containerList.add(atomContainer2);

            SpinSystem spinSystem = new MagneticPropTool().calculateHSpinSystem(atomContainer);
            SpinSystem spinSystemH = new MagneticPropTool().calculateAverageHSpinSystem(containerList);
            //this is just a basic test
            // TODO think of a better test
            for (int i = 0; i < spinSystemH.getNbOfChemEquivalentSpinGroups(); i++){
                Assert.assertTrue("Something changed in the Bolztmann averaging",
                        (spinSystem.getMagneticTensors()[i]-spinSystemH.getMagneticTensors()[i])<0.2);
            }


        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }


}
