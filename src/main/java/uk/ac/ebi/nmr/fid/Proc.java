package uk.ac.ebi.nmr.fid;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 14/01/2013
 * Time: 14:02
 * To change this template use File | Settings | File Templates.
 */
public class Proc {

    private static int windowFunctionType;      //wdw               window function type
    private static int phasingType;             //ph_mod            phasing type
    private static int f1DetectionMode;         //mc2               F1 detection mode
    private static double zeroFrequency;        //sf                frequency of 0 ppm
    private static double lineBroadning;        //lb                line broadening (in Hz?)
    private static double gbFactor;             //gb                GB-factor
    private static double zeroOrderPhase;       //phc0              zero order phase
    private static double firstOrderPhase;      //phc1              first order phase
    private static double ssb;                  //ssb               sine bell shift
    private static double ssbSine;              //ssbSine           sine bell shift
    private static double ssbSineSquared;       //ssbSineSquared    sine bell shift

    // obtained after reading the FID but could it be obtained from the Bruker?
    private static int transformSize;           //si                transform size (complex)    bruker_read::get_fid
    private static double dwellTime;            //dw                dwell time (in s)           bruker_read::get_fid
    private static double hertzPerPoint;        //hzperpt                                       bruker_read::get_fid
    private static double ppmPerPoint;          //ppmperpt                                      bruker_read::get_fid
    private static double spectraWidthHertz;    //sw_h                                          bruker_read::get_fid
    private static double offset;               //offset                                        bruker_read::get_fid

    //TODO consider moving this to the Fourier Transformed class or processing class....
    // variables required later...
    private static int tdEffective;             //td_eff        apodization::transform::do_fft
    private static int leftShift=0;             //leftshift     ft_settings_dialog::ft_settings_dialog
    private static int shift;                   //j             apodization::transform::do_fft
    private static int increment;               //i             apodization::transform::do_fft

    public Proc(Acqu acquisition) throws Exception {

        // set the size for the fourier Transform
        // perhaps I should check first if I can use the data from the Proc file..??
        if (acquisition.getAquiredPoints() < 1*1024) this.transformSize =1024;
        else if (acquisition.getAquiredPoints()<= 2 * 1024) this.transformSize=2*1024;
        else if (acquisition.getAquiredPoints() <= 4 * 1024) this.transformSize=4*1024;
        else if (acquisition.getAquiredPoints() <= 8 * 1024) this.transformSize=8*1024;
        else if (acquisition.getAquiredPoints() <= 16 * 1024) this.transformSize=16*1024;
        else if (acquisition.getAquiredPoints() <= 32 * 1024) this.transformSize=32*1024;
        else if (acquisition.getAquiredPoints() <= 64 * 1024) this.transformSize=64*1024;
        else if (acquisition.getAquiredPoints() <= 128* 1024) this.transformSize=128*1024;
        else if (acquisition.getAquiredPoints() <= 256* 1024) this.transformSize=256*1024;
        else this.transformSize=512 * 1024;
        //set the dwell time (in s) to display the timeline of the fid (dw is distance between points)
        if(acquisition.getSpectralWidth() == 0 | acquisition.getTransmiterFreq() == 0)
            throw new Exception ("Some acquisition parameters are null");
        this.dwellTime = 1.0/(2 * acquisition.getSpectralWidth() * acquisition.getTransmiterFreq());
        this.hertzPerPoint = acquisition.getSpectralWidth() * acquisition.getTransmiterFreq() / transformSize;
        this.ppmPerPoint = acquisition.getSpectralWidth() / transformSize;
        this.spectraWidthHertz = acquisition.getSpectralWidth() * acquisition.getTransmiterFreq();
        this.offset = (acquisition.getTransmiterFreq()- zeroFrequency) * 1.0e06 +
                (acquisition.getSpectralWidth()* acquisition.getTransmiterFreq()) / 2.0;

        // set the position where the shift starts???
        switch (acquisition.getAcquisitionMode()) {
            case DISP:
            case SIMULTANIOUS:
                shift=2*leftShift;
                break;
            case SEQUENTIAL:
                shift=leftShift;
                break;
            default:
                break;
        }

        //set the number of acquired points we are going to work with
        tdEffective=(acquisition.getAquiredPoints()<=2*transformSize)?
                acquisition.getAquiredPoints():                     // normal case
                2*transformSize;                                  // fid data is truncated in the nonsense case

    }

    public static int getTransformSize() {
        return transformSize;
    }

    public static void setTransformSize(int transformSize) {
        Proc.transformSize = transformSize;
    }

    public static int getWindowFunctionType() {
        return windowFunctionType;
    }

    public static void setWindowFunctionType(int windowFunctionType) {
        Proc.windowFunctionType = windowFunctionType;
    }

    public static int getPhasingType() {
        return phasingType;
    }

    public static void setPhasingType(int phasingType) {
        Proc.phasingType = phasingType;
    }

    public static int getF1DetectionMode() {
        return f1DetectionMode;
    }

    public static void setF1DetectionMode(int f1DetectionMode) {
        Proc.f1DetectionMode = f1DetectionMode;
    }

    public static double getZeroFrequency() {
        return zeroFrequency;
    }

    public static void setZeroFrequency(double zeroFrequency) {
        Proc.zeroFrequency = zeroFrequency;
    }

    public static double getLineBroadning() {
        return lineBroadning;
    }

    public static void setLineBroadning(double lineBroadning) {
        Proc.lineBroadning = lineBroadning;
    }

    public static double getGbFactor() {
        return gbFactor;
    }

    public static void setGbFactor(double gbFactor) {
        Proc.gbFactor = gbFactor;
    }

    public static double getZeroOrderPhase() {
        return zeroOrderPhase;
    }

    public static void setZeroOrderPhase(double zeroOrderPhase) {
        Proc.zeroOrderPhase = zeroOrderPhase;
    }

    public static double getFirstOrderPhase() {
        return firstOrderPhase;
    }

    public static void setFirstOrderPhase(double firstOrderPhase) {
        Proc.firstOrderPhase = firstOrderPhase;
    }

    public static double getSsb() {
        return ssb;
    }

    public static void setSsb(double ssb) {
        Proc.ssb = ssb;
        // if ssb is given in degrees this converts it to the inverse of coefficient in front to the Pi
        // 360 = 2 Pi => angle/180 = coefficient (e.g. 360/180 = 2 Pi)
        // This variables are used in the apodizationTool method sine and Co.
        if (ssb >= 1)                //convert Bruker convention to degrees
            ssb = 180.0 / ssb;
        else
              ssb = 0.0;
        // I do not get this... (see source code)
        Proc.ssbSine = ssb;
        Proc.ssbSineSquared = ssb;
    }

    public static double getSsbSine() {
        return ssbSine;
    }

    public static double getSsbSineSquared() {
        return ssbSineSquared;
    }

    public static int getTdEffective() {
        return tdEffective;
    }

    public static void setTdEffective(int tdEffective) {
        Proc.tdEffective = tdEffective;
    }

    public static int getLeftShift() {
        return leftShift;
    }

    public static void setLeftShift(int leftShift) {
        Proc.leftShift = leftShift;
    }

    public static int getShift() {
        return shift;
    }

    public static void setShift(int shift) {
        Proc.shift = shift;
    }

    public static int getIncrement() {
        return increment;
    }

    public static void setIncrement(int increment) {
        Proc.increment = increment;
    }

    public static double getDwellTime() {
        return dwellTime;
    }
}
