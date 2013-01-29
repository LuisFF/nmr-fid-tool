package uk.ac.ebi.nmr.fid;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 14/01/2013
 * Time: 14:01
 * To change this template use File | Settings | File Templates.
 */
public class Acqu {


    //TODO use an enum for parameters that have a limited set of options such as aquisition mode
    private static double transmiterFreq;               //sfo1
    private static double decoupler1Freq;               //sfo2
    private static double decoupler2Feq;                //sfo3
    private static double freqOffset;                   //o1 (Hz)
    private static double spectralWidth;                //sw            (ppm)
    private static int aquiredPoints;                   //td            acquired points (real + imaginary)
    private static int dspDecimation;                   //decim         DSP decimation factor
    private static int dspFirmware;                     //dspfvs        DSP firmware version
    private static int dspGroupDelay;                   //grpdly        DSP group delay
    private static int byteOrder;                       //bytorda       byte order
    private static int filterType;                      //digmod        filter type
    private static int numberOfScans;                   //ns            number of scans
    private static boolean integerType;                 //dtypa         data type (0 -> 32 bit int, 1 -> 64 bit double)
    private static String pulseProgram;                 //pulprog       pulse program
    private static String observedNucleus;              //nuc1          observed nucleus
    private static String instrumentName;               //instrum       instrument name
    private static String solvent;                      //solvent       solvent
    private static String probehead;                    //probehead     probehead
    private static String origin;                       //origin        origin
    private static String owner;                        //owner         owner
    private static AcquisitionMode acquisitionMode;      //aq_mod        acquisition mode
    private static FidData fidType;                     //fid_type      define in class data_par



    public Acqu() {

    }

    public static double getTransmiterFreq() {
        return transmiterFreq;
    }

    public static void setTransmiterFreq(double transmiterFreq) {
        Acqu.transmiterFreq = transmiterFreq;
    }

    public static double getDecoupler1Freq() {
        return decoupler1Freq;
    }

    public static void setDecoupler1Freq(double decoupler1Freq) {
        Acqu.decoupler1Freq = decoupler1Freq;
    }

    public static double getDecoupler2Feq() {
        return decoupler2Feq;
    }

    public static void setDecoupler2Feq(double decoupler2Feq) {
        Acqu.decoupler2Feq = decoupler2Feq;
    }

    public static double getFreqOffset() {
        return freqOffset;
    }

    public static void setFreqOffset(double freqOffset) {
        Acqu.freqOffset = freqOffset;
    }

    public static double getSpectralWidth() {
        return spectralWidth;
    }

    public static void setSpectralWidth(double spectralWidth) {
        Acqu.spectralWidth = spectralWidth;
    }

    public static int getAquiredPoints() {
        return aquiredPoints;
    }

    public static void setAquiredPoints(int aquiredPoints) {
        Acqu.aquiredPoints = aquiredPoints;
    }

    public static int getDspDecimation() {
        return dspDecimation;
    }

    public static void setDspDecimation(int dspDecimation) {
        Acqu.dspDecimation = dspDecimation;
    }

    public static int getDspFirmware() {
        return dspFirmware;
    }

    public static void setDspFirmware(int dspFirmware) {
        Acqu.dspFirmware = dspFirmware;
    }

    public static int getDspGroupDelay() {
        return dspGroupDelay;
    }

    public static void setDspGroupDelay(int dspGroupDelay) {
        Acqu.dspGroupDelay = dspGroupDelay;
    }

    public static int getByteOrder() {
        return byteOrder;
    }

    public static void setByteOrder(int byteOrder) {
        Acqu.byteOrder = byteOrder;
    }

    public static AcquisitionMode getAcquisitionMode() {
        return acquisitionMode;
    }

    public static void setAcquisitionMode(int acquisitionMode) {
        for (AcquisitionMode mode : AcquisitionMode.values())
            if (mode.type == acquisitionMode)
                Acqu.acquisitionMode = mode;
    }
    public static void setAcquisitionMode(AcquisitionMode mode) {
        Acqu.acquisitionMode=mode;

    }

    public static int getFilterType() {
        return filterType;
    }

    public static void setFilterType(int filterType) {
        Acqu.filterType = filterType;
    }

    public static int getNumberOfScans() {
        return numberOfScans;
    }

    public static void setNumberOfScans(int numberOfScans) {
        Acqu.numberOfScans = numberOfScans;
    }

    public static String getPulseProgram() {
        return pulseProgram;
    }

    public static void setPulseProgram(String pulseProgram) {
        Acqu.pulseProgram = pulseProgram;
    }

    public static String getObservedNucleus() {
        return observedNucleus;
    }

    public static void setObservedNucleus(String observedNucleus) {
        Acqu.observedNucleus = observedNucleus;
    }

    public static String getInstrumentName() {
        return instrumentName;
    }

    public static void setInstrumentName(String instrumentName) {
        Acqu.instrumentName = instrumentName;
    }

    public static boolean is32Bit() {
        return integerType;
    }

    public static void set32Bit(boolean is32Bit) {
        fidType=(is32Bit)?FidData.INT32:FidData.DOUBLE;
        Acqu.integerType = is32Bit;
    }

    public static String getSolvent() {
        return solvent;
    }

    public static void setSolvent(String solvent) {
        Acqu.solvent = solvent;
    }

    public static String getProbehead() {
        return probehead;
    }

    public static void setProbehead(String probehead) {
        Acqu.probehead = probehead;
    }

    public static String getOrigin() {
        return origin;
    }

    public static void setOrigin(String origin) {
        Acqu.origin = origin;
    }

    public static String getOwner() {
        return owner;
    }

    public static void setOwner(String owner) {
        Acqu.owner = owner;
    }

    public static FidData getFidType() {
        return fidType;
    }

    public enum AcquisitionMode {
        SEQUENTIAL      (1),
        SIMULTANIOUS    (2),
        DISP            (3);

        private final int type;

        private AcquisitionMode(int type){
            this.type=type;
        }

        private double getType(){return type;};
    }

    public enum FidData {
        INT32       (1),
        DOUBLE      (2),
        FLOAT       (3),
        INT16       (4);

        private final int type;

        private FidData(int type){
            this.type=type;
        }

        private double getType(){return type;};
    }
}
