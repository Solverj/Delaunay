import java.util.*;
import java.util.concurrent.*;
import easyIO.*;


// Delaunay Triangulation by Arne Maus 2014
// reimplementation testing speed and different algorithms

class DelaunayAlg1 {
    static DelaunayAlg1 d;
	public static void main(String [] args) {

		if (args.length != 4 ) {
			System.out.println(" use: >java DelaunayAlg1 <n num points> <points  per box>  <debugfil> <dataOutput>");
		} else {
			fil = args[2];
			data = new Out(args[3]);
			 d =new DelaunayAlg1();
			d.PerformAlg(args);

		}
	} // end main

    DT dt;
    static String fil ;
    static Out ut, data ;

	static void println(String s) {
		ut= new Out(fil,true);
		ut.outln(s) ;
		ut.close();
		System.out.println(s);
	}

	static void print(String s) {
			ut= new Out(fil,true);
			ut.out(s) ;
		    ut.close();
			System.out.print(s);
	}

	void outputDT() {
			// output n datapoints sorted)
				data.outln(n+"");
				// output CoHull
				data.out("Cohull, num: "+ theCoHull.numElem+"points:");
				for (int i = 0; i < theCoHull.numElem; i++){
						data.outln(theCoHull.dtEdges[i]+" ");
				 }

				// output Datapoints
					for (int i = 0;i < n; i++){
						data.outln(""+i + " " +x[i] +" "+y[i]);
				   }

				 // output Delaunay neghbours to al points - first is nearest neighbour
				   int [] nb;
				 	for (int i = 0;i < n; i++){
						nb = delEdges[i];
						data.out(""+i+ " num= "+nb.length+": ");

						for (int j =0; j < nb.length;j++){
							data.out(nb[j]+" ");
						}
						data.outln();

					}


			   data.close();
    }// end outputPoints


	void PerformAlg(String [] args){
		// do Delaunay triangulation ----------------------------------------(start main prog, perfor alg)---------------

			n = Integer.parseInt(args[0]);
			NUM_PER_BOX = Integer.parseInt(args[1]);

		    // println("\nDelaunay triangulering med: " + n + " punkter og ca."+NUM_PER_BOX+" punkter per boks. "+
		        //        ", fil:"+args[2]+"\n");

		    long t;
		    double cTime,dtTime,complTime;
		    int numRest=0;

		    //1)  algorithm
		    t = System.nanoTime();
		    initiate();
		    cTime = (double) (System.nanoTime()-t)/1000000.0;
		   // println("\nTid: Initiate:" + Format.align(cTime,10,3)+" msek, dvs. "+Format.align((cTime*1.0/n),10,5)+ " msek. per punkt");
			 //    dumpPoints("Initiate");



		   // 2) find convex hull
			t = System.nanoTime();
			int numCholl =  cohull(allNB);
			cTime = (double) (System.nanoTime()-t)/1000000.0;
			println("\nTid: Cohull fant:"+theCoHull.numElem+" punkter på cohull. Tid:" + Format.align(cTime,10,3)+" msek, dvs. "
			  +Format.align((cTime*1.0/n),10,6)+ " msek. per punkt");
			//if (n < DebugLimit)			dumpPoints("Etter cohull");

		  //  dt = new DT (d);

		 // 3) Find near neighbours and dt edges - do not need the convex hull
		    t = System.nanoTime();
		    numNNeigh= doAlg1Delaunay(theCoHull);
		    dtTime= (double) (System.nanoTime()-t)/1000000.0;
		    // se min Bit-artikkel
		    if (n < DebugLimit)dumpPoints("Etter Alg1 - alle ");

		    // 4a) Check consistency - matching neighbours
		    //sjekkDK (true);
		    int missing = sjekkDK(true);


		    // 4) some statistics
		    antIndrePunkter = n - theCoHull.numElem;
		    totAntKanter = 2*theCoHull.numElem + 3*antIndrePunkter -3;
		    totAntTriangler=  theCoHull.numElem + 2*antIndrePunkter -2;
		    println("Ant pkt:"+n+", antIndrePkt:"+antIndrePunkter+", totAntKanter:"+totAntKanter+", totAntTriangler:"+totAntTriangler);

		    println("\nDT-kanter funnet:"+ (numNNeigh/2)+ " unike med " +NUM_PER_BOX+" gj.snitt punkter per box på "  +
		                        Format.align(dtTime,10,3)+"msek , \ndvs. "+Format.align(dtTime/n,8,4)+
		                        " millisek. per punkt, og:" +Format.align(2*(numNNeigh*1.0/(2*n)),10,5)+
		                       " snitt kanter. per punkt.\n Prosent av alle funnet:"+
		                        Format.align((100.0*(numNNeigh/2)/totAntKanter),10,4)+"%"+
		                        "\n prosent funnet (inkl missing):"+
		                        Format.align((100.0*((numNNeigh/2)+missing/2)/totAntKanter),10,4)+"%");

	  	   for (int i = 0; i < n; i++) {
			   if (delEdges [i].length < minDK ){minDK = delEdges [i].length; minPDK =i;}
			   if (delEdges [i].length > maxDK ){maxDK =delEdges [i].length; maxPDK =i;}
		   }
		   print("Min # DE :"+minDK+" for point:"+minPDK);
		   println(", max # DE :"+maxDK+" for point:"+maxPDK);
		   println("Cocirkulartitet (4 eller flere punkter) funnet:"+ Format.align((numCoCircular/4.0),10,2)+"  ganger.");

	//	   dumpPoints("Etter nærmeste nabo");


		 if (n < 2000)    dt = new DT (d);

    }//----------------------------------------------------------------( end PerformAlg, start Data structure) --------

	 /** Data Structure - Delaunay model: */
	 int n ;           // number of points
	 int numNNeigh =0; // number of neigbours tested avg
	 final static int  SCALING = 9,
	                 MIN_NUM_FOUND = 5;  // upper average limit
	 int NUM_PER_BOX,
	     MAX_X,
		 MAX_Y,
		 NUM_BOX,        // number of boxes (in x and y - direction)
		 BOX_SIDE,
		 NUM_X_BOXES,
		 NUM_Y_BOXES,
		 numCoCircular =0,
		 searchSize;     // sourroung box when searching for close neighbours (= 3,5,7,..)
	int DebugLimit = 60 ;
     double EPSILON;        // Tolerance for comparing double value (cocirculat case)
	 int [] x, 				// x-values of points
			y; 				// y-value of points
	 int [][]box;			// grid structure in both x- and y over area. Approx NUM_PER_BOX in ech box;
	 int [][]delEdges ;     //delEdges[i][] = All Delaunay Edges to point 'i'
	 //int [] numDE;         // numDE[i] =  number of Delaunay Edges to point i
	 int minDK = 60, minPDK =0, maxDK =0, maxPDK;  // min and max DK found by nearesr neighbours

	 int antIndrePunkter;         // = n - theCoHull.numElem;
	 int totAntKanter;            // = 2*theCoHull.numElem + 3*antIndrePunkter -3;
	 int totAntTriangler;         // =  theCoHull.numElem + 2*antIndrePunkter -2;
	 int   numK4 =0;                 // number of DK funnet av Th. 4

	 HashMap <Integer,Integer> coHullMap;
	 NList theCoHull = new NList();           // The points of the convex hull
//	       theCoCircular = new NList();       // A list of the points that are coCircular (initiall not ordered)
     // end data model DT of n points ---------------------------------------------(end data model)---------------------
    int [] allNB = new int [100] ;   // assumed way beyond max number og neigbours to a node


	 // Debugging tool
     void dumpPoints(String s) {
		 println ("Dump av DataPoints: "+s);
		 for (int i=0;i < n; i++){
		   print(" Point i:"+i+" ("+x[i]+","+y[i]+") med DK til:");
		   if (delEdges[i] != null ) {
			   for (int j =0; j < delEdges[i].length; j++){
				 if(j >0) print(",");
				 print(delEdges[i][j]+"");
			   }
		   }
		   println("");
		 }
	  }// end dumpPoints

	  // Sjekk konistens av DK-kanter (Hvis a-b, så b-a
	  // ikke sjekk på kryssende kanter
	  int sjekkDK (boolean output) {
		  int numErr =0, numPrinted =0, maxPrinted =10;
		  for (int p = 0; p < n; p++) {
			  // for each delaunay nerghbour pn to p, check that p is Del-neighb to pn
			  for (int pn : delEdges[p]) {
				  boolean ok = false;
				  for (int i : delEdges[pn]){
					  if (i == p) { ok = true; break;}
				  }
			      if (! ok) {
					  numErr ++;
					  if (output && numPrinted < maxPrinted ){
						    numPrinted++;
						    println ("\n"+numPrinted+")DT-error, punkt:"+p+"("+x[p]+","+y[p]+") har :" + pn +"("+x[pn]+","+y[pn]+
						            ") som DelaunayKant, men ikke motsatt");
						    print("\n  D-neighb for "+p+": ");
							for (int j =0; j < delEdges[p].length; j++){
								 if(j >0) print(",");
								 print(delEdges[p][j]+"");
							}
							print("\n  D-neighb for "+pn+": ");
							for (int j =0; j < delEdges[pn].length; j++){
								 if(j >0) print(",");
								 print(delEdges[pn][j]+"");
							}
							print("\n");
							if ( numPrinted ==  maxPrinted)
							println ("...............etc.................");
						} // erroroutput
				  } // error
			   } // end all D neighbours to p
		   } // end p - all points
           println("Number of not matching Delaunay neghbours:"+numErr);
		   return numErr;
	  }// end sjekkDK


	// debugging
	 void dumpBox(String s){
		 print("\nBOX-dump:"+s+", BOX_SIDE:"+BOX_SIDE+", NUM_X_BOXES:"+NUM_X_BOXES+", MAX_X:"+MAX_X);
		 for (int xx = 0; xx <= NUM_X_BOXES; xx++){
			 print("\nbox["+xx+"][]:");
			 for (int yy = 0; yy<= NUM_Y_BOXES;yy++){
				 print(box[xx][yy]+",");
			 }
		 }
		 println("");
	 }// end dumpBox

    // Initiate data structure ------
	void initiate(){
	// ------------------------------------------------------------------(start initialize data structure)---------------

	   MAX_X = (int) Math.sqrt(n)*SCALING +1; // must have room for non-data points in integer x,y-plane
	   MAX_Y = MAX_X;
	   EPSILON = 10.0/MAX_X;
	   int numP=0;
	   long mb1,mb2;

	   Random r = new Random(1237);
	   x = new int [n];
	   y = new int [n];
	   Runtime runtime = Runtime.getRuntime();
	   delEdges = new int[n][];

	   NUM_BOX = n/NUM_PER_BOX+1;
	   // make NUM_BOX a square
	   while (! (((int) Math.sqrt(NUM_BOX))*((int) Math.sqrt(NUM_BOX)) == NUM_BOX))
	         NUM_BOX++;
	   NUM_X_BOXES = (int)Math.sqrt(NUM_BOX);
	   NUM_Y_BOXES= NUM_X_BOXES;
	   BOX_SIDE = (int) (MAX_X/NUM_X_BOXES);
	   if (NUM_X_BOXES*BOX_SIDE < MAX_X) BOX_SIDE++;
	   box = new int[NUM_X_BOXES+1][NUM_Y_BOXES+1];  // need dummy last elements

	   // make BOX_SIDE at least as big as BOX_SIDE*NUM_BOX >= MAX_X
	   while ( BOX_SIDE*NUM_X_BOXES < MAX_X) BOX_SIDE++;

	  // cohull = new boolean[n]; // default 'false'
	  coHullMap = new HashMap<Integer,Integer>((int) Math.sqrt(n)+10); // apprx size og coHull

	  // 1)  initiate with random numbers -
	  // might here alternativly read datapoints from file
	  new NPartallsPunkter(n).fyllArrayer(x,y);

println("** NUM_BOX:"+NUM_BOX+",BOX_SIDE:" + BOX_SIDE+", MAX_X/Y:" + MAX_X+ ", NUM_X/Y_BOXES:"+NUM_X_BOXES+ ", MIN_NUM_FOUND:"+MIN_NUM_FOUND);

	  // 2) Sort on first x and then for each x-stripe on y
	  //    create boxes and their borders
	  radix2(x,y,0,x.length); // sort on x and move corresponding y-values

	 // dumpPoints("After first x-sort");

	  int low =0, high  =0;
	  int xVal = BOX_SIDE;
	  int boxIndex=0;
	  box[0][0] =0 ; // start of first box
	//  box[NUM_X_BOXES] [NUM_Y_BOXES+1]= n; // terminate last box

	  int [] bbox;   // dummy for box[i]; i.e those points with similar x-values

	  // find the boxes
      for (int i = 0; i <= NUM_X_BOXES; i++) {   // test i < NUM_X_BOXES if extra row and collumn in box is removed
		    bbox = box[i];
		    boxIndex =0;


			// find new x-stripe
			low = high;
			while( high < n && x[high] < xVal) high++;

//			println("X-stripe:"+i+", low:"+low+", high:"+high);

			//sort that x-stripe on y
			radix2(y,x,low,high);

			// find all boxes in that stripe
			int  yVal = BOX_SIDE;

			bbox[boxIndex] = low;  // found start of  box[boxIndex]

			int yInd = low;
			for (int j = 0; j < NUM_Y_BOXES; j++) {
			   while (yInd < high && y[yInd] < yVal) {
		//		  println(" +++  P : "+yInd+"("+x[yInd]+","+y[yInd]+")");
				   yInd++;
			   }
      // println("------------");

			   bbox[++boxIndex] = yInd;
			   yVal += BOX_SIDE;

			} // end new y-box  (for j)
			xVal += BOX_SIDE;

		    bbox[NUM_Y_BOXES] = yInd;
//			    println("------------");
		} // end find boxes


		if (n < DebugLimit) dumpPoints("After sort");
		if (n < DebugLimit) dumpBox("After Sort");
	} // -----------------------------------------------------------------------( end initiate, start cohull) ----------;


	  /*******************************************************************************
	   *  find the convex hull of all poits, Add the edges on the cohull to the set
	   *   of Delaunay Edges and return number of points on 'cohull'.
	   *******************************************************************************/
	   int cohull(int []allNB){
		   int minx =0,maxx=0,miny=0,maxy=0;
		   for (int i = 1; i<n; i++){
			   if (x[i] < x[minx]) minx =i;
			   if (x[i] > x[maxx]) maxx =i;
			   if (y[i] < y[miny]) miny =i;
			   if (y[i] > y[maxy]) maxy =i;
		   }
		   theCoHull.put(minx);
		   coHullMap.put(minx,minx);

		// println(" minx:"+minx+", miny:"+miny+", maxx:"+maxx+", maxy:"+maxy);
		   if (miny != minx){
			   theCoHull.putAfter(miny,minx);
			   coHullMap.put(miny,miny);
			   cohullRec(minx,miny);
		    } // new point miny
		   if (maxx != miny) {
			   theCoHull.putAfter(maxx,miny);
			   coHullMap.put(maxx,maxx);
			   cohullRec(miny,maxx);
	       } // new point maxx
		   if (maxy != maxx && maxy != minx){
			   theCoHull.putAfter(maxy,maxx);
			   	coHullMap.put(maxy,maxy);
			   cohullRec(maxx,maxy);
		   }
		   if (maxy != minx) cohullRec(maxy,minx);
           else cohullRec(maxx,minx);

		   return theCoHull.numElem;
	   }// end cohull

		/************************************************************************************
		 * Find  a possible next point minP between p1-p2 closest to p1 on the convex hull;
		 * recurse on p1-minP  and  minP-p2 .
		 *************************************************************************************/
		void cohullRec(int p1,int p2) {
	 // println("\n\n++ cohullrec param, p1;"+p1+",p2:"+p2);
			   // finding line ax+by+c =0 from p1 to p2
			   int [] linep1p2 = new int[3];
			   linep1p2 = line (linep1p2,p1,p2);
			   int a = linep1p2[0], b = linep1p2[1], c = linep1p2[2], p2Save = p2;
			   int minD2=0;  // max squared distance (negative) for point outside p1-p2
			   int d2;
			   int minP=-1;
			   int numNeg =0;
			   // open boxes
			   int numx = Math.abs(x[p1]-x[p2])/BOX_SIDE +1;
			   // for all relevant x-stripes


//if (n < DebugLimit)
// println("\n   - cohullRec tester linje : p1 "+p1+"("+x[p1]+","+y[p1]+") - p2:"
//		            +p2+"("+x[p2]+","+y[p2]+   "): "+a+"x + "+b+"y + "+c +"=0");

			   int lowx, highx, lowy, highy;
			   if ( x[p1] <= x[p2]) {
				   //from left to right - bottom
				   lowx = x[p1]/BOX_SIDE;
				   highx= x[p2]/BOX_SIDE;
				   lowy = 0;
				   highy = Math.max(y[p1]/BOX_SIDE,y[p2]/BOX_SIDE);
			   }else{
				   // from right to left - top side
				   lowx = x[p2]/BOX_SIDE;
				   highx= x[p1]/BOX_SIDE;
				   lowy= Math.min(y[p1]/BOX_SIDE,y[p2]/BOX_SIDE);
				   highy = NUM_Y_BOXES-1;
			   }


	//if (n < DebugLimit)
	// println("     --cohull lowx:"+lowx+", highx:"+highx+", lowy:"+lowy+", highy:"+highy);

			   for ( int xx =lowx ; xx <= highx; xx++){
				     // For every relevant y-box on these stripes
					  for (int i= box[xx][lowy]; i < box[xx][highy+1]; i++){

						  if ( insideRect(p1,p2,i)) {
							  // for all other points in this box
							  d2 = a*x[i] +b*y[i] +c;
			  // println("     Testing p '"+i+"':"+x[i]+","+y[i]+", verdi:"+d2);
			                  if (d2 <= 0) numNeg++;
                              if (  d2 <= minD2 ){
								  // found better  point i
								  minP = i;
								  minD2 = d2;
								 // println(" *** Found better point on cohull p:"+minP);
							  } // end better point

						  }// end test point 'i' in box(xx,yy)
					  }// end open box(xx,yy)
			   } // end x-stripe

			   if (numNeg > 0 ) {
				  // found new point outside OR ON line p1-p2

				  theCoHull.putAfter(minP,p1);
				  coHullMap.put(minP,minP);     // add to HashMap
                  if ( numNeg > 1) {
					 cohullRec (p1,minP);
				     cohullRec(minP,p2);
				 }

			   }
			}// end ------------cohullRec----------------------

			//** returns true if point i ifinside rectangle defined by p1 and p2 */
			boolean insideRect(int p1,int p2, int i) {
				boolean xOK = false, yOK=false;

				// println ( "        insideRect( p1:"+ p1+"("+x[p1]+","+y[p1]+")"+", p2:"+
				 //                                p2 +"("+x[p2]+","+y[p2]+")"+", i:"
				 //                                +i+"("+x[i]+","+y[i]+")"+");)");

				if (x[p1] < x[p2]) xOK = x[p1] <= x[i] && x[i] <= x[p2];
				else xOK = x[p1] >= x[i] && x[i] >= x[p2];

				if (y[p1] < y[p2]) yOK = y[p1] <= y[i] && y[i] <= y[p2];
				else yOK = y[p1] >= y[i] && y[i] >= y[p2];

				return xOK & yOK & (p1 !=i) & (p2!=i);
			}// end inside p1p2Rectangle

// ----------------------------------------------------------------------------------- (end finding cohull) --------------

     /*      /** translates from user x- and y-ccordinates to box-adress /
		   int boxAdress(int xx, int yy) {
			   // returns box adress of point p (xx,yy)
			   return  (( BOX_SIDE*NUM_Y_BOXES )/xx + yy/BOX_SIDE);
		   } // end boxAdress
       */
		 /** translates from box row (x) , collum (y) number to array-indexes in x[] and y[] */
		   int boxAdress2(int stripex, int rowy) {
			   // returns box adress of row, collumn
			   return stripex*NUM_Y_BOXES + rowy;
		   } // end boxAdress

           /** translates from user (double) x- and y-ccordinates to array-indexes in x[] and y[] */
		   int boxAdress3(double xx, double  yy) {
			   // returns box adress of point p (xx,yy)
			   return (int) (( BOX_SIDE*NUM_Y_BOXES )/xx + yy/BOX_SIDE);
		   } // end boxAdress




 // ----------------------------------------------------------(start findClosestNeighbourEdges)---------------------
       double [] line2 (double [] line, double ax, double ay, double bx, double by) {
			// returnes the line equation cx+dy+e =0: line[0]x+ line[1]y+ line[2]=0,
			// line from  a to b
			//double [] line = new double[3];
			line[0] = ay - by;
		    line[1] = bx - ax;
	        line[2] = by*ax-ay*bx;
			return line;
		}// end line

	   // 2.jan 2015
	   int [] line (int [] line, int ax,int ay, int bx, int by) {
				// returnes the line equation cx+dy+e=0 : line[0]x+ line[1]y+ line[2]=0,
				// line from  a to b
			//	int [] line = new int[3];
				line[0] = ay - by;
			    line[1] = bx - ax;
		        line[2] = by*ax-ay*bx;
				return line;
		}// end line

	  // line from point a to point b
	  int [] line(int [] line, int a,int b) {
				// returnes the line equation cx+dy+e=0 : line[0]x+ line[1]y+ line[2]=0,
				// line from  pint a to b
		//		int [] line = new int[3];
				line[0] = y[a] - y[b];
				line[1] = x[b] - x[a];
				line[2] = y[b]*x[a]-y[a]*x[b];
				return line;
		}// end line

		int distToLine(int [] line, int point){
			return line[0]*x[point] +line[1]*y[point] + line[2];
		} // end distToLine

		double distToLine(int [] line, double xx, double yy){
				return line[0] *xx +line[1]*yy + line[2];
		} // end distToLine

        double [] midNormal (int a, int b) {
			 // returns the midnormal from line A - B on the form ret[0] + ret[1]y+ ret[1]
             double []ret = new double[3];

			 // (mx,my) midpoint an a-b,
			 int mx = (x[a]+x[b])/2;  // OK, x and y are even numbers
			 int my = (y[a]+y[b])/2;
			 double p3x ,p3y = y[a];
			 if ( x[b] == x[a]){ p3x = 2*mx;p3y=my;}  // nmidnorm is parallel with x-axis
			 else if (y[b] == y[a]) {p3x =mx; p3y=2*my;} // midnorm parallel y-axis
			 else p3x =(x[b]*x[b] - x[a]*x[a] + Math.pow(y[b] - y[a],2))/(2.0*(x[b]-x[a]));
			// return line (mx,my,p3x,y[a]);
            ret = line2 (ret,mx,my,p3x,p3y);
          // println("--midNormal on: "+a+"-"+b+" , with LineEq "+ Format.align(ret[0],6,1)
            //         +"*x +"+Format.align(ret[1],6,1)+" *y + " + Format.align(ret[2],6,1));
			return ret;
		 }// end midNormal

		 double [] midNormCrossing(double [] midAB, double [] midAC, double [] crossing){
			 // returns the crossing, x=crossing[0], y = crossing [1] of two lines midAB and midAC)
			 double div = midAC[0]*midAB[1]-midAB[0]*midAC[1];
			 // if (Math.abs(div) < 0.00001) return false;  // Paralell lines will not accur
			 crossing[0] = (midAC[1]*midAB[2]-midAC[2]*midAB[1])/div; // x-value
			 crossing[1] = (midAC[2]*midAB[0]-midAC[0]*midAB[2])/div; // y -value
			 return crossing;
		 }// end midNormCrossing

	      boolean cAbovepB (int c, int [] linepB) {
			  // test which side of a line pB point c is (>0 to the left = 'above'
		      return (linepB[0]*x[c] + linepB[1]*y[c] + linepB[2]) > 0;
	      }// end findCfromBase

	   long dist2 (int x, int y){
		   return x*x+y*y;
	   } // end dist2 int

	   double dist2 (double x, double y){
	   		   return x*x+y*y;
	   } // end dist2 double


		boolean pointJinCirclePI (int mix, int miy, int j, int rip2){
			// Has point j shorter radius to circle p-i than the squared radius of th circle
			return dist2(mix-x[j],miy-y[j]) < rip2;
		}// end pointJinCirclePI


      /******************************************************************
       *	Alg 1: Finn først nærmeste nabo b for punkt a, så
       *           finn neste nabo C til A, dvs. trtekanten ABC ved:
       *		1. Først først den C med minst (lavest, gjerne  neg verdi) på
       *		   skjæring mellom midtnormalene AB og AC blant alle punkter
       *			inne i de bokser som dekker A og B.
       *		2. Test så om funnet C er inne i sirkel som er inkludert i
       *		   de boksene det er søkt i  og 'over' AB.
       *		   a) Ok , C er funnet
       *		   b) NEI, utvid søkebokser til de dekker Sirkel ABC
       * 				og søk i resten av dem (etter evt bedre C)
       ********************************************************************/
       int doAlg1Delaunay (NList chull) {
		   int b,stop, num =0;

		   // a) first find all neighbours to cohull
		   for (int i = 0; i <= chull.numElem-1; i++){
			   if ( i == 0)  stop  =  chull.dtEdges[ chull.numElem -1]; else stop  = chull.dtEdges[i-1];
			   if (i == chull.numElem-1) b = chull.dtEdges[0]; else b= chull.dtEdges[i+1];
			   num += findRestOfNeighbours(chull.dtEdges[i], b, stop, allNB);

		   }

		   // b) then find first closest neighbour, and then rest of neighbours to all points a
		   for (int a = 0; a < n; a++) {
		      if ( !coHullMap.containsKey(a)) {
				  b = closestNeigbourTo(a);
				  allNB[0] = b;
	//if (n < DebugLimit)println("NNeighb to:"+a+" is:"+b);
				  num += findRestOfNeighbours(a,b,b, allNB); // stop when finding b again
			  }
		   }
		   return num;
	   }// end doAlg1Delaunay



       // find and add all  neighbours to ''a'from start DK ab. Return when next neighbor is stop
       // When calling this method, only b has been added as neighbour to point a.
       int  findRestOfNeighbours (int a, int b, int stop, int[] allNB) {
		   int numNB =0;
		   int c =b;
		   boolean cohullCase = (stop!=b);

		   // println("\n\nfindRestOfNeighbours from :"+a+" to "+ b+", stop:"+stop+"------------");

		   // find all neighbour  anti clockwise round a
		   do{
			  b =c;
			  allNB[numNB++] = c;
			  c = findNextNeighbour (a,b);
		   } while (c != stop);

		   if (cohullCase ) { // cohull case
		     numNB++;
		   }

		   // set in neighbourlist
		   int[] nbToA  = new int[numNB];

		   if (cohullCase ) { // cohull case
			 allNB[numNB-1] = stop;
		   }

		   for (int i = 0; i < numNB; i++){
			   nbToA [i] = allNB[i];
		   }
		   delEdges[a] = nbToA;
           return numNB;
	   } // end findRestOfNeighbours;

       // ************************************************************************
	   // find next Delaunay neighbor c from line a-b - ie. next D-triangle abc
	   //***************************************************************************
	   int findNextNeighbour(int a, int b) {
		   int c = -1; // not legal value
		   double cr2 = 0;     // Dummy start , squared dist c to S
		   int mx = (x[a]+x[b])/2, my = (y[a]+y[b])/2;

		   // find searccircle S center (sx,sy) along midtnormal a-b
		   double dx = ((y[b] - y[a])*5.0)/17, dy = ((x[b]-x[a])*5.0)/17; // angle c appr = 60 deg.
		    int []lineAB = new int[3];
		    lineAB = line(lineAB,x[a],y[a],x[b],y[b]);

		//  if (n < DebugLimit) println("\n*FindNextNeighbour, a:"+a+" b :"+b +", Equation:" +
		//        Format.align(lineAB[0],6)+"* x +"+Format.align(lineAB[1],6)+"* y +"+
		 //        Format.align(lineAB[2],10));

		   // bestem midtnormal a-b
		   double [] midNormAB = midNormal(a,b),
		             midNormAP;
           double sx = mx, sy = my;

           // Central loop - loop until at least one potential c (p) is found - pick best
		   do {
			   sx = sx-dx;
			   sy = sy+dy;
		       double sr2 = dist2(x[a]-sx, y[a]-sy); // squared radius for S
		       int numBox=0;
		       // find search block area that includes S

		       // r = minste avstand for S til box-vegg
			   double xMin = Math.min(sx%BOX_SIDE, BOX_SIDE-(sx%BOX_SIDE)),
					  yMin = Math.min(sy%BOX_SIDE, BOX_SIDE-(sy%BOX_SIDE)),
					  min = Math.min(xMin,yMin),
					  max = (BOX_SIDE*numBox)+min,
				      maxR2 = max*max;
			   double minDistToAB = 10000000.0;
			   double []crossing = new double[2];
			   int [] coCirc = new int[10];
			   int coCircular =0;

		       while (maxR2 < sr2) {
				   max += BOX_SIDE;
				   maxR2= max*max;
				   numBox++;
			   }

//println(" - A1 , sx:"+Format.align(sx,8,2)+", sy:"+Format.align(sy,8,1)+", max:"+Format.align(max,6,2)+
//                  ", maxR2:"+Format.align(maxR2,8,2)+", sr2:"+Format.align(sr2,8,1)+ ",\n BOX_SIDE:"+ BOX_SIDE+", numBox:"+numBox );


		       // finn lowx, lowy, highx, highy;
			   int lowx =  Math.max(0,((int)sx/BOX_SIDE)-(numBox+1));
			   int highx = Math.min(((int)sx/BOX_SIDE)+(numBox)+1, NUM_X_BOXES-1);
			   int lowy =  Math.max(0,(int)sy/BOX_SIDE -(numBox+1));
			   int highy = Math.min((int)sy/BOX_SIDE +(numBox+1), NUM_Y_BOXES-1);

// println(" -- A2 , lowx:"+lowx+", highx:"+highx+", lowy:"+lowy+", highy:"+highy );


		       // løp gjennom omådet og finn den beste c som
		       for (int xx = lowx; xx <= highx; xx++) {

				   for (int p = box[xx][lowy]; p < box[xx][highy+1];  p++) {
			//	 println(" --- B1, xx:"+xx+", p:"+p);
					   if (p != a && p!=b) {
						   double radp2 = (x[p]-sx)*(x[p]-sx)+(y[p]-sy)*(y[p]-sy);

			//			   println(" ----C0 , Test om p:" +p+"  med r2:"+radp2+", i sirkel med  maxR2:"+maxR2);

						   if ( radp2 <= sr2 &&  cAbovepB(p, lineAB)) {
								// point p within search circle &&  'above' line a-b
								midNormAP = midNormal(a,p); // Find midnormal A-C

									 //println(" ------ C1, Fant brukbar p:"+p+", nåværende c:"+c);

								crossing = midNormCrossing (midNormAB,midNormAP,crossing);
								// normal, all points above ab are not colinear points with ab
								double dist = distToLine(lineAB, crossing[0], crossing[1]);

								//	 println(" ---- B2, p:"+p+", radp2:"+Format.align(radp2,10,1)+", dist:"
								//	          +Format.align(dist,20,1)+ ", minDistToAB:"+Format.align(minDistToAB,10,1));

								if (Math.abs(dist - minDistToAB) < EPSILON  && c != p){
									  //  cocircular points c and p
									  if(coCircular == 0) {
										// must add a,b and best c found so far (for later calculations)
										coCirc[coCircular++] = c;
									  }
									           // println(" ******* CoCirc, FANT CoCirc p:"+p+" to c:"+c+
									            //       ", diff dist:"+ Math.abs(dist - minDistToAB));
								    coCirc [coCircular++] = p;
								} else if ( dist < minDistToAB) {
									  // found new better c
									  c = p;
													  // println(" --------- C2, FANT OK c:"+c);
									  minDistToAB = dist;
									  coCircular = 0;
								} // end not cocircular
							 }// end usable point p
						}// end p!= a,b
					} // end new candidate point p for c in search circle
               } // end xx

               // must find right CoCircular point if needed
				if (coCircular > 0) {
				  // must choose the same cocirc point as all orher A,b on that circle
				  c = findCoCirc(a,b,coCircular, coCirc);
	//println(" ++++++++++++Fant cocirc trekant, a:"+a+", b:"+b+", c:"+c);
	              numCoCircular++;
				  coCircular = 0;
				}
			   dx +=dx;
			   dy+= dy;
		   } while ( c < 0);

		   return c;

	   }// end findNextNeighbour

	   // find c in coCircular case, Algorithm:
	   // let c1 = closest point to origo, and if two points have
	   //     same distance, the c1 = those of these two with smallest x-value
	   // if (c1 != a && c != b) return c1 else
	   //  if (c1 == a  ) take as c the next point  to left of a; else
	   // if (c1 == b) take as c the nest point to the to right of b
	   //    Sub-alg finding the one point p to the left of a (similar to the right of b):
	   // try all lines a-p , p one of the other cocirc points nort a,b,
	   // and the  p1 without without any point p =1 a,b above a-p1 is then c
	   //
	   int findCoCirc(int a, int b,int coCircular, int []coCirc){
		   long dMin2 = MAX_X*MAX_X, d2;
		   int [] lineC1p=new int[3];
		   int c1=-1,p,q, num=0, i,j;

		   // add a & b for first test
		   coCirc[coCircular] =a;
		   coCirc[coCircular+1] = b;

		   //find p closest to origo
		   for ( i = 0; i < coCircular+2; i++) {
			   p = coCirc[i];
			   d2 = x[p]*x[p];
			   if (c1 < 0 ||(d2 == dMin2 && y[p] < y[c1])|| d2 < dMin2){
				   c1 = p; dMin2 =d2;
			   }
		   } // found candidate c1 closest to origo

		   if (c1 == b) {
				// find as c the next point to the left of a
			    for ( i = 0; i < coCircular; i++) {
					p = coCirc[i];
					lineC1p = line(lineC1p,c1,p);
					num=0;
					for ( j = 0; j < coCircular; j++) {
							q = coCirc[j];
							if (p!=q){
	//	 println("--++findCoCirc A  a:"+a+",b:"+b+", c1 funnet:"+c1+", og prøver p:"+p+", og q:"+q);
								if (  cAbovepB (q, lineC1p)) {
								  // point i har another point j above it - p is not our point
                      //             println("      A "+q+" above "+c1+"-"+p);
								   break ;  // test next point i
								}
							} // test all points q != p
					 } // end j, test all cocirc points q for above line: a-p
					    if (j == coCircular) {
							// reached end of j-loop without fing a point above c1-p
					        c1 = p; // p is our point because no points q above it
					        break; // this i-loop is finished
						}
			    } // end i
			} else if ( c1== a ) {
				// find as c the next point to the right of b
				for ( i = 0; i < coCircular; i++) {
					p = coCirc[i];
					lineC1p = line(lineC1p,c1,p);
						for ( j = 0; j < coCircular; j++) {
							q = coCirc[j];
							if (p!=q){
		// println("--++findCoCirc B  a:"+a+",b:"+b+", c1 funnet:"+c1+", og prøver p:"+p+", og q:"+q);
								if ( ! cAbovepB (q, lineC1p)) {
								  // point i har another point j below it - p is not our point
						//		   println("    B   "+q+" NOT above "+c1+"-"+p);
								   break ;  // test next point i
							   } // end if
						   }//  end test all points q != p
						} // end j,  test all cocirc points  for above line: a-p
						if (j == coCircular) {
							c1 = p; // p is our point because no points q below it
							break; // this i-loop is finished
						}
				 } // end i
			 } // end if c1==a

			 return c1; //
 	   } // end findCoCirc



	   int closestNeigbourTo(int p){
		   // start with b =0 a legal value
		   int b=0;
		   int xp = x[p], yp =y[p];
           long bDist2 = MAX_X*MAX_X;
           int numBox =0;
           long maxR2 = 0;

		   int xMin = Math.min(xp%BOX_SIDE, BOX_SIDE-xp%BOX_SIDE),
			   yMin = Math.min(yp%BOX_SIDE, BOX_SIDE-yp%BOX_SIDE),
			   min = Math.min(xMin,yMin),
			   max;

           do  {
			   numBox +=1;
			   // while not found closest
			   int lowx =  Math.max(0,((int)xp/BOX_SIDE)-(numBox+1));
			   int highx = Math.min(((int)xp/BOX_SIDE)+(numBox)+1, NUM_X_BOXES-1);
			   int lowy =  Math.max(0,(int)yp/BOX_SIDE -(numBox+1));
			   int highy = Math.min((int)yp/BOX_SIDE +(numBox+1), NUM_Y_BOXES-1);

			   max = (BOX_SIDE*numBox)+ min;
			   maxR2 = max*max;
//println(" \n NNa closestNeigbourTo;"+p+", lowx:"+lowx+", highx:"+highx+", lowy:"+lowy+
//     ", highy:"+highy+ "numBox:"+numBox +", p:"+p);

			   for (int xx = lowx; xx <= highx; xx++) {
				   for (int j = box[xx][lowy]; j < box[xx][highy+1];j++) {
					   if (j!= p){
						 //// println("   --- testing j:"+j);
						   long d2 = dist2(x[j]-xp,y[j]-yp);
						   if ( d2 < bDist2) {
							   // found new, better b
							   b = j;
							    // println("   ---   FOUND j:"+j +", bDist2:"+bDist2);
							   bDist2 = d2;
						   }
					   }// end p != j
				   } // end j
			   } // end i

		   }  while ( bDist2 > maxR2) ;// test if b is within circle
//println(" NNb closestNeigbourTo;"+p+"is: "+b);

		   return b;
	   } // end closestNeigbourTo

// end finding DelEdges by nearest neighbours ------------------------------------(end findClosestNeighbourEdges)-----------



       // --------------------------------------------------------------------(start NList = liste for cohull og coCircular)----
	   class NList {
		 static final int START_NUM_NEIGHBOURS = 7;
		 int numElem =0;
		 int p ;  // this list is owned by point p
		 int [] dtEdges;

		 NList() {
			 dtEdges = new int[START_NUM_NEIGHBOURS];
		 }

		 void dump(String s) {
			 println ("Dump av liste:"+s);
			 for (int i=0;i < numElem; i++){
			   int t = dtEdges[i];
			  println(" Point t:"+t+" ("+x[t]+","+y[t]+")");
			 }
		 }// end dump

		 /** put newDE after suc */


		 void  putAfter(int newDE, int suc){
//			 println(" * putAfter nytt pkt:"+  newDE+" after:"+suc);
			 if (numElem == dtEdges.length) {
				 // expand dtEdges
				 int[] b = new int [(int)(dtEdges.length *2)]; // double each iteration
				 for (int i = 0; i < dtEdges.length; i++){
					 b[i] = dtEdges[i];
				 }
				 dtEdges = b;
			  } // end expand array
			  // all cases:
			  int i = numElem-1;
			  // make room for new member
			  while (i >=0 && dtEdges[i] != suc){
				  dtEdges[i+1]= dtEdges[i];
				  i--;
			  }
			  dtEdges[i+1] = newDE;
			  numElem++;
//		  println(" --- A Cohull:puttet nytt punkt nr. "+newDE+" etter:"+suc+", numKanter.:" + numElem);
		 } // end put new Delaunay edge after 'suc'

		/** put number of point p in in nList[i] at end - then  p-i is a Delaunay Edge */
		void  put(int newDE){
//			println(" * put nytt pkt:"+  newDE);
			 if (numElem == dtEdges.length) {
				 // expand dtEdges
				 int[] b = new int [(int)(dtEdges.length *1.41)]; // double every second iteration
				 for (int i = 0; i < dtEdges.length; i++){
					 b[i] = dtEdges[i];
				 }
				 dtEdges = b;
			  } // end expand array
			  // all cases:
			   dtEdges[numElem++] = newDE;

	      //println(" --- put: nytt punkt nr. "+newDE+" inn i NList");
		 } // end put new Delaunay edge

		 void delete(int p) {};

		/** sortEdgesConterClockwise */
		void sortEdgesConterClockwise() {
		};
	 } // ----------------------------------------------------------------end NList , start sorting--------------


     /** SEKVENSIELL    insertsort (int[],int[])
	  * Sort on a[low..high] and move y[low..high]  same way */
     static void  insertSort(int [] a, int []y, int low, int high) {
		     int i, t,ty;

		      for (int k = low ; k < high; k++){
		          t =  a[k+1];
		          ty = y[k+1];
		          i = k;

		          while (i >= low && a[i] > t) {
		             a[i+1] = a[i];
		             y[i+1] = y[i];
		             i--;
		           }
		          a[i+1] = t;
		          y[i+1] = ty;
		     }
    } // end insertSort -int


     /** SEKVENSIELL    insertsort (double[],int[])
	  * Sort on a[low..high] and move y[low..high]  same way */
	  static void insertSort(double a[],int[] y,int low, int high) {
	     int i, ty;
	     double t;

	  //   println("A insertSort (double) low:"+low+", high:"+high);

	      for (int k = low ; k < high; k++) {
	          t =  a[k+1];
	          ty = y[k+1];
	          i = k;
		//	 println(" - B insertSort (double) k:"+k+", i:"+i+", a[i]:"+a[i]+", t:"+t+", ty:"+ty);

	          while (i >= low && a[i] > t) {
		//		  println("   - C insertSort (double) k:"+k+", i:"+i+", a[i]:"+a[i]+", t:"+t+", ty:"+ty);

	             a[i+1] = a[i];
	             y[i+1] = y[i];
	             i--;
	           }
	          a[i+1] = t;
	          y[i+1] = ty;
	     } // end k
	 }// end double,int insertSort

    /** SEKVENSIELL    insertsort (long, int)
	  * Sort on a[low..high] and move y[low..high]  same way */
     static void  insertSort(long [] a, int []y, int low, int high) {
		     int i,ty;
		     long t;

		     for (int k = low ; k < high; k++){
		          t =  a[k+1];
		          ty = y[k+1];
		          i = k;

		          while (i >= low && a[i] > t) {
		             a[i+1] = a[i];
		             y[i+1] = y[i];
		             i--;
		           }
		          a[i+1] = t;
		          y[i+1] = ty;
		     }//end k

    } // end insertSort -long


    /** SEKVENSIELL    R A D I X  ***********************
     * Sort on a[low..high> and move y[low..high>  same way */
	  static void radix2(int [] a, int []y, int low, int high) {
	//	  println(" Sort radix2 with low:"+low+", high:"+high+ ", a.length:" + a.length);
		  if (low < high) {
			 if ( high - low < 100 ) {
				 insertSort(a, y, low, high-1) ;
		     } else {

				  // 2 digit radixSort on a[low..high>, let y do same moves
				  int max = a[low], numBit = 2;

				  for (int i = low+1 ; i < high ; i++)
					   if (a[i] > max) max = a[i];

				  while (max >= (1<<numBit) )numBit++;

				  int bit1 = numBit/2,
					  bit2 = numBit-bit1;

				  int[] b = new int [high-low+1];
				  int[] yb = new int[high-low+1];

				  radixSort( a,b,low,high,y,yb,bit1, 0);
				  radixSort( b,a,low,high,yb,y,bit2, bit1);
			  } // end radix-sort
		  } // end insertsort (or no sort)
	  } // end radix2

		 /** Sort a[] on one digit ; number of bits = maskLen, shiftet up
			 ‘shift’ bits */
		 static void radixSort ( int [] a, int [] b,int low, int high, int []y, int yb[], int maskLen, int shift){
			  int  acumVal = 0, j, n = a.length;
			  int mask = (1<<maskLen) -1;


			  int [] count = new int [mask+1];
			  int ind;
			  // new variables because of no backcopy of b to a
			  int alow;
			  if (a.length > b.length) {
				  // normal first  digit, soritng orig 'a' to shorter 'b'
				  alow =0;
			  }else {// second digit, soritng shorter 'b' to orig 'a'
				  alow =low;
			  }

			 // a) count=the frequency of each radix value in a
			  for (int i = low; i < high; i++)
				 count[(a[i-alow]>> shift) & mask]++;

			 // b) Add up in 'count' - accumulated values
			  for (int i = 0; i <= mask; i++) {
				   j = count[i];
					count[i] = acumVal;
					acumVal += j;
			  }

			 // c) move numbers in sorted order a to b
			  for (int i = low; i < high; i++) {
				 ind =  count[(a[i-alow]>>shift) & mask]++;
				 b[ind+alow]  = a[i-alow];
				 yb[ind+alow] = y[i-alow];
			  }// end move x (in a[]) and y

       }// end Sekvensiell radixSort **************************************

}// end class Delaunay

