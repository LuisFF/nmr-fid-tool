package uk.ac.ebi.nmr;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Chemical equivalence groups are defined based on molecular symmetry.
 * The molecular symmetry can be obtained from the HuLuIndexTool, i.e. atoms with the same index are symmetrical
 *
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 10/08/2012
 * Time: 14:31
 * To change this template use File | Settings | File Templates.
 */
public class SpinSystem {

    private IAtomContainer atomContainer;
    private int nbOfChemEquivalentSpinGroups;
    private double energy;

    private int [] numberOfSpins;
    private double [] magneticTensors;
    private double [][] constantJMatrix;

    public SpinSystem(){

    }

    public SpinSystem(int numberOfClasses) {
        numberOfSpins = new int[numberOfClasses];
        magneticTensors = new double[numberOfClasses];
        constantJMatrix = new double[numberOfClasses][numberOfClasses];
        nbOfChemEquivalentSpinGroups = numberOfClasses;
    }

    public int getNbOfChemEquivalentSpinGroups() {
        return nbOfChemEquivalentSpinGroups;
    }

    public void setNbOfChemEquivalentSpinGroups(int nbOfChemEquivalentSpinGroups) {
        this.nbOfChemEquivalentSpinGroups = nbOfChemEquivalentSpinGroups;
    }

    public double[] getMagneticTensors() {
        return magneticTensors;
    }

    public void setMagneticTensors(double[] magneticTensors) {
        this.magneticTensors = magneticTensors;
    }

    public double[][] getConstantJMatrix() {
        return constantJMatrix;
    }

    public void setConstantJMatrix(double[][] constantJMatrix) {
        this.constantJMatrix = constantJMatrix;
    }

    public int[] getNumberOfSpins() {
        return numberOfSpins;
    }

    public void setNumberOfSpins(int[] numberOfSpins) {
        this.numberOfSpins = numberOfSpins;
    }

    public double getEnergy() {
        return energy;
    }

    public void setEnergy(double energy) {
        this.energy = energy;
    }

    public IAtomContainer getAtomContainer() {
        return atomContainer;
    }

    public void setAtomContainer(IAtomContainer atomContainer) {
        this.atomContainer = atomContainer;
    }
}
