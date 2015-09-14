import java.util.Random;
/**
* Class NPunkter for aa finne n tilfeldige, ulike  punkter i x,y-planet
* ver 7.mai 2015
************************************************************************/
public class NPartallsPunkter {
	Random r;
	int n;
	byte [] bitArr;
	static int maxXY, xShift =3;
	int scaleFactor = 3;  // scaleFactor * scaleFactor * n= antall mulige punkter i planet (her: 4*n)
	final  int [] bitMask ={1,2,4,8,16,32,64,128};

	NPartallsPunkter(int n) {
		this.n =n;
		maxXY = Math.max(10,(int) Math.sqrt(n) * scaleFactor); // største X og Y verdi
		while ((1<<xShift) < maxXY) xShift++;
		xShift = xShift - 3;                    // 8 bits per byte
		bitArr = new byte[(maxXY<<xShift |(maxXY>>3))  + 5];
		r = new Random(123);
	}

	private void setUsed(int x, int y) {
      bitArr[(x<<xShift) | (y >>3)] |= bitMask[(y&7)];
	}

	private boolean used (int x, int y) {
		return  (bitArr[(x<<xShift) | (y >>3)] & bitMask[y&7]) != 0;
	}

	public void fyllArrayer(int [] x, int[] y) {
		int next =0;
		int xval, yval,maxHalve = maxXY/2;;
		while (next < n) {
			do{
				xval = r.nextInt(maxHalve)+1;
				 xval <<=1;          // now an even number
				yval = r.nextInt(maxHalve)+1;
				yval <<=1;  // now an even number
			} while (used (xval, yval));
			x[next] = xval;
			y[next] = yval;
			setUsed (xval,yval);
			next++;
		} // next point
	}// end fyllArrayer

} // end class NPunkter
