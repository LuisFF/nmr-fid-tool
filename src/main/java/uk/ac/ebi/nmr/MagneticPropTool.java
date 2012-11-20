package uk.ac.ebi.nmr;

import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.config.AtomTypeFactory;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.invariant.EquivalentClassPartitioner;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import uk.ac.ebi.nmr.io.G09OutputReader;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 28/08/2012
 * Time: 10:28
 * To change this template use File | Settings | File Templates.
 */
public class MagneticPropTool {

    public static final String EQUIVALENT_CLASS = "org.openscience.cdk.nmr.magneticproptool:equivalent.class";

    private final static double BOLTZMANN_CONSTANT = 1.3806488E-23; // J K^-1
    private final static double BOLTZMANN_CONSTANT2 = 8.6173324E-5; // eV K^-1
    private final static double HARTREE_ENERGY = 4.35974434E-18; // J
    private final static double HARTREE_ENERGY2 = 27.21138505; // eV

//    private IAtomContainer originalAtomContainer = null;
//    private IAtomContainer atomContainerWithoutH = null;
//    private HashMap<Integer, List<IAtom>> partitionMap = null;
//    private HashMap<IAtom, Integer> atomClassMap = null;
//    private HashMap<IAtom, IAtom> mapContainers = null;
//    private SpinSystem spinSystemH = null;


    public MagneticPropTool() {


    }

    /**
     * map the atoms between the container without hydrogens and the original container.
     * The objects are no longer the same and one way to map them back is through the 3D coordinates.
     */
    private HashMap<IAtom, IAtom> mapAtoms(IAtomContainer atomContainerWithoutH,
                                           IAtomContainer originalAtomContainer) {
        HashMap<IAtom, IAtom> mapContainers = new HashMap<IAtom, IAtom>();
        //To change body of created methods use File | Settings | File Templates.
        for (IAtom atom1 : atomContainerWithoutH.atoms()) {
            for (IAtom atom : originalAtomContainer.atoms()) {
                if (atom.getPoint3d().equals(atom1.getPoint3d())) {
                    mapContainers.put(atom1, atom);
                    break;
                }
            }
        }
        return mapContainers;

    }

    /**
     * perform the atom partitioning into equivalent classes in order to define the topological equivalent spins.
     */
    private HashMap<Integer, List<IAtom>> partitionAtoms(IAtomContainer atomContainerWithoutH, HashMap<IAtom, IAtom> mapContainers) {
        HashMap<Integer, List<IAtom>> partitionMap = new HashMap<Integer, List<IAtom>>();


        int[] partitions = null;
        try {
            // partition the atoms into equivalent classes
            partitions = new EquivalentClassPartitioner(atomContainerWithoutH)
                    .getTopoEquivClassbyHuXu(atomContainerWithoutH);

            /*
            reset the number of partitions in order to consider only carbon atoms.
            Map the old partitionMap value to a new one
            do the mapping already with the original IAtom objects
             */
            HashMap<Integer, Integer> mapNewPartitions = new HashMap<Integer, Integer>();
            partitions[0] = 1;

            for (int i = 1; i < partitions.length; i++) {
                if (atomContainerWithoutH.getAtom(i - 1).getSymbol().equals("C")) {
                    if (!mapNewPartitions.containsKey(partitions[i])) {
                        mapNewPartitions.put(partitions[i], partitions[0]);
                        partitions[0]++;
                        partitionMap.put(mapNewPartitions.get(partitions[i]), new ArrayList<IAtom>());
                    }
                    mapContainers.get(atomContainerWithoutH.getAtom(i - 1))
                            .setProperty(EQUIVALENT_CLASS, mapNewPartitions.get(partitions[i]));
                    partitionMap.get(mapNewPartitions.get(partitions[i]))
                            .add(mapContainers.get(atomContainerWithoutH.getAtom(i - 1)));
                }

//                System.out.println(partitions[i]+ " " +
//                        atomContainerWithoutH.getAtom(i-1).getImplicitHydrogenCount()+ " "+
//                        atomContainerWithoutH.getAtom(i-1).getAtomTypeName());
            }

//        SMSDNormalizer.convertExplicitToImplicitHydrogens(atomContainerWithoutH);
//        this.atomContainerWithoutH = atomContainer;
        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return partitionMap;
    }

    /**
     * calculates the average spin system for protons assuming a Boltzmann distribution based on
     * the molecular energy in Hartree-Fock approximation.
     *
     * @param atomContainers
     * @return
     * @throws Exception
     */
    public SpinSystem calculateAverageHSpinSystem(List<IAtomContainer> atomContainers) throws Exception {
        //TODO abstract the boltzmann weighted averaging

        List<SpinSystem> spinSystemList = new ArrayList<SpinSystem>();
        String inchi = null;
        // create the spinSystems
        for (IAtomContainer container : atomContainers) {
            if (inchi == null)
                inchi = InChIGeneratorFactory.getInstance().getInChIGenerator(container).getInchi();
            // check if we are talking about the same molecule
            if (!(InChIGeneratorFactory.getInstance().getInChIGenerator(container).getInchi().equals(inchi)))
                System.err.println("The AtomContainer does not represent the same molecule\n"+
                        InChIGeneratorFactory.getInstance().getInChIGenerator(container).getInchi()+"\n"+
                inchi);
            spinSystemList.add(calculateHSpinSystem(container));
        }

        // average the magnetic tensor and and coupling constants of the spinsystem list;
        // I'll assume that the atom order did not change...
        double[] averageShift = null;
        double[][] averageCoupling = null;
        double energySum = 0;
        SpinSystem averageSpinSystem = new SpinSystem();
        double minEnergy = 0;
        for (SpinSystem spinSystem : spinSystemList) {
            minEnergy = (spinSystem.getEnergy() < minEnergy) ? spinSystem.getEnergy() : minEnergy;
        }

        for (SpinSystem spinSystem : spinSystemList) {
            // declare arrays and pass invariant properties
            if (averageShift == null) {
                // this things should be invariant.... cross-fingers...
                averageSpinSystem.setNbOfChemEquivalentSpinGroups(spinSystem.getNbOfChemEquivalentSpinGroups());
                averageSpinSystem.setNumberOfSpins(spinSystem.getNumberOfSpins());
                averageSpinSystem.setAtomContainer(spinSystem.getAtomContainer());
                System.out.println(spinSystem.getNbOfChemEquivalentSpinGroups());
                // declare averageShift and averageCoupling arrays
                averageShift = new double[spinSystem.getNbOfChemEquivalentSpinGroups()];
                averageCoupling = new double[spinSystem.getNbOfChemEquivalentSpinGroups()]
                        [spinSystem.getNbOfChemEquivalentSpinGroups()];
            }
            // calculate the energy dependent factor specific to the conformation of the spin system
            // temperature set to 300K
            System.out.println("Relative energy of the conformer: " + (spinSystem.getEnergy() - minEnergy));
            System.out.println("Exponent term in Bolztamnn distribution: "
                    + (-(spinSystem.getEnergy() - minEnergy) * HARTREE_ENERGY / (BOLTZMANN_CONSTANT * 300)));
            double factor = Math.exp(-((spinSystem.getEnergy() - minEnergy) * HARTREE_ENERGY / (BOLTZMANN_CONSTANT * 300)));
            // sum up the energy dependent factor
            energySum += factor;
            // sum the magnetic tensors (and coupling constants) weighted by the energy dependent factor
            for (int i = 0; i < spinSystem.getMagneticTensors().length; i++) {
                if (spinSystem.getNumberOfSpins()[i] > 0) {
                    averageShift[i] += spinSystem.getMagneticTensors()[i] * factor;
                    for (int j = i + 1; j < spinSystem.getMagneticTensors().length; j++) {
                        if (spinSystem.getConstantJMatrix()[i][j] != 0) {
                            averageCoupling[i][j] += spinSystem.getConstantJMatrix()[i][j] * factor;
                            averageCoupling[j][i] = averageCoupling[i][j];
                        }
                    }
                }
            }
        }

        // perform the averaging step using the sum of the energy dependent factor
        for (int i = 0; i < averageShift.length; i++) {
            if (averageShift[i] > 0) {
                averageShift[i] = averageShift[i] / energySum;
                for (int j = i + 1; j < averageShift.length; j++) {
                    averageCoupling[i][j] = averageCoupling[i][j] / energySum;
                    averageCoupling[j][i] = averageCoupling[i][j];
                }
            }
        }
        // assign the average magnetic tensor and coupling constant matrix to the object variables
        averageSpinSystem.setMagneticTensors(averageShift);
        averageSpinSystem.setConstantJMatrix(averageCoupling);

        return averageSpinSystem;
    }

    /**
     * Calculates the average tensor and coupling constants for each of the protons in the same spin system.
     * The spin systems is defined based on topological equivalence.
     *
     * @param atomContainer
     * @return
     * @throws Exception
     */
    public SpinSystem calculateHSpinSystem(IAtomContainer atomContainer) throws Exception {
        // first thing to do is remove the explicit hydrogens
        // the atoms can be retrieved by comparing their 3d coordinates
        // yes, at this point the atoms should have 3d coordinates
        IAtomContainer originalAtomContainer = atomContainer;
        IAtomContainer atomContainerWithoutH = SilentChemObjectBuilder.getInstance()
                .newInstance(IAtomContainer.class, atomContainer.clone());
        SpinSystem spinSystemH = null;
        try {
            AtomTypeFactory factory = AtomTypeFactory.getInstance("org/openscience/cdk/config/data/jmol_atomtypes.txt"
                    , atomContainerWithoutH.getBuilder());
            for (IAtom atom : atomContainerWithoutH.atoms())
                factory.configure(atom);

            CDKHueckelAromaticityDetector.detectAromaticity(atomContainerWithoutH);
            SmilesGenerator smiler = new SmilesGenerator(true);
            System.out.println("original compound: " + smiler.createSMILES(atomContainerWithoutH));
            System.out.println("counts before explicit: " +
                    AtomContainerManipulator.getTotalHydrogenCount(atomContainerWithoutH) + " " +
                    atomContainerWithoutH.getAtom(7).getAtomTypeName() + " " +
                    atomContainerWithoutH.getAtom(1).getImplicitHydrogenCount() + " " +
                    atomContainerWithoutH.getAtom(1).getBondOrderSum());

            // remove the hydrogens so that the EquivalenteClassPartitioner works properly
            atomContainerWithoutH = AtomContainerManipulator.removeHydrogensPreserveMultiplyBonded(atomContainerWithoutH);
            //        AtomContainerManipulator.convertImplicitToExplicitHydrogens(atomContainerWithoutH);
            System.out.println("compound without hydrogens: " + smiler.createSMILES(atomContainerWithoutH));
            System.out.println("counts after explicit hydrogens: " +
                    AtomContainerManipulator.getTotalHydrogenCount(atomContainerWithoutH) + " " +
                    atomContainerWithoutH.getAtom(1).getAtomTypeName() + " " +
                    atomContainerWithoutH.getAtom(1).getImplicitHydrogenCount() + " " +
                    atomContainerWithoutH.getAtom(1).getBondOrderSum());
            // maps the atoms of the atomContainerWithoutH to the originalAtomContainer using the 3D coordinates
            HashMap<IAtom, IAtom> mapContainers = mapAtoms(atomContainerWithoutH, originalAtomContainer);
            // perform that atom partitioning into equivalent classes
            // not that the list of atoms contains already the atoms from the originalAtomContainer
            HashMap<Integer, List<IAtom>> partitionMap = partitionAtoms(atomContainerWithoutH, mapContainers);
            // first check how many carbons have hydrogens and store their indices
            List<Integer> indicesCWithH = new ArrayList<Integer>();
            System.out.println("Partitioning number: "+ partitionMap.size());
            spinSystemH = new SpinSystem(partitionMap.size());
            int[] numberOfSpins = new int[partitionMap.size()];
            double[] avgMagneticTensor = new double[partitionMap.size()];
            double[][] avgCouplingConstant = new double[partitionMap.size()][partitionMap.size()];

            for (int i = 0; i < partitionMap.size(); i++) {
                for (IAtom atom : partitionMap.get(i + 1)) {
                    // iterate through the hydrogens attached to the carbon and update the magnetic tensor value
                    // assuming that they belong to the same equivalent class because they are attached to the same carbon
                    for (IAtom spins : originalAtomContainer.getConnectedAtomsList(atom)) {
                        if (spins.getSymbol().equals("H")) {
                            numberOfSpins[i]++;
                            avgMagneticTensor[i] += (Double) spins.getProperty(G09OutputReader.MAGNETIC_TENSOR);
                            // iterate through all hydrogens that are coupled with this one and check to which class their
                            // carbon belongs to and update the coupling value to that equivalent classes
                            HashMap<IAtom, Float> coupledAtoms = (HashMap<IAtom, Float>) spins
                                    .getProperty(G09OutputReader.COUPLING_CONSTANTS);
                            for (IAtom coupleAtom : coupledAtoms.keySet()) {
                                if (coupleAtom.getSymbol().equals("H")) {
                                    // check the number of coupled atoms
                                    if (originalAtomContainer.getConnectedAtomsCount(coupleAtom) > 1)
                                        throw new Exception("Hydrogen is connected to more than one atom");
                                    // check the number of coupled atoms
                                    if (originalAtomContainer.getConnectedAtomsCount(coupleAtom) == 0)
                                        throw new Exception("Hydrogen is disconected: " +
                                                (new SmilesGenerator().createSMILES(originalAtomContainer)));
                                    // get the class of the hydrogen and add the coupling constant is not from the same class
                                    IAtom connectedAtom = originalAtomContainer.getConnectedAtomsList(coupleAtom).get(0);
                                    if (connectedAtom.getSymbol().equals("C") && (Integer) connectedAtom.getProperty(EQUIVALENT_CLASS) != (i + 1))
                                        avgCouplingConstant[i][((Integer) connectedAtom.getProperty(EQUIVALENT_CLASS)) - 1] +=
                                                coupledAtoms.get(coupleAtom).doubleValue();
                                }
                            }
                        }
                    }

                }

            }
            /*
            average the magnetic tensors
            make ConstantJ matrix squared
            */
            for (int i = 0; i < numberOfSpins.length; i++) {
                if (numberOfSpins[i] > 0)
                    avgMagneticTensor[i] = avgMagneticTensor[i] / numberOfSpins[i];
                if (i > 0) {
                    for (int j = 0; j < i; j++) {
                        if (numberOfSpins[j] > 0 && numberOfSpins[i] > 0) {
                            avgCouplingConstant[i][j] = avgCouplingConstant[i][j] / (numberOfSpins[j] * numberOfSpins[i]);
                            avgCouplingConstant[j][i] = avgCouplingConstant[i][j];
                        }
                    }
                }
            }

            /*
           Average the magnetic tensor and coupling constant values
           Make coupling matrix symmetric
            */
            spinSystemH.setAtomContainer(originalAtomContainer);
            spinSystemH.setNumberOfSpins(numberOfSpins);
            spinSystemH.setConstantJMatrix(avgCouplingConstant);
            spinSystemH.setMagneticTensors(avgMagneticTensor);
            spinSystemH.setEnergy((Double) originalAtomContainer.getProperty(G09OutputReader.STRUCTURE_ENERGY));
            System.out.println("Energy of the spin system: " + spinSystemH.getEnergy());
        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return spinSystemH;
    }


    public void calculateEquivClassShift4C() {

        // TODO implement


    }


}
