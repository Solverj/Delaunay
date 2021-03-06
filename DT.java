import javax.swing.*;
import java.awt.*;
import java.awt.geom.Line2D;

class DT extends JFrame{
	DelaunayAlg1 d;

	DT(DelaunayAlg1 d){
		this.d =d;
		size =600;
		margin =50;
	    scale =size*1.0/d.MAX_X;
		  setTitle("DelaunayDemo, num points:"+d.n);
		  grafen = new Graph();
		 // add (grafen, BorderLayout.CENTER);

		  getContentPane().add(grafen, BorderLayout.CENTER);
		  setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		  pack();
		  setVisible(true);

		// angir foretrukket størrelse på dette lerretet.
		setPreferredSize(new Dimension(d.MAX_X+2*margin,d.MAX_Y+2*margin));
	}


	Graph grafen;
	int size , margin;
	double scale ;


    class Graph extends JPanel{
	    int SIZE ;

	    void drawPoint(int p, Graphics g) {
                 SIZE = 7;

			    //  g.drawString(""+p,xDraw(d.x[p]-SUB),yDraw(d.y[p]+SUB));
			    if (d.n < d.DebugLimit) g.drawString(p+"("+d.x[p]+","+d.y[p]+","+d.z[p]+")",
			          xDraw(d.x[p],0),yDraw(d.y[p],SIZE/2+1));
			    else if (d.n < d.DebugLimit*2) g.drawString(p+"",xDraw(d.x[p],0),yDraw(d.y[p],SIZE/2+1));
				 g.drawOval (xDraw(d.x[p],SIZE/2),yDraw(d.y[p],SIZE/2),SIZE,SIZE);
				 g.fillOval (xDraw(d.x[p],SIZE/2),yDraw(d.y[p],SIZE/2),SIZE,SIZE);
	     }

		 Graph() {
			 setPreferredSize(new Dimension(size+2*margin+10,size+2*margin+10));
		 }

	     int   numX = (int) Math.ceil(d.MAX_X/d.BOX_SIDE),
		       numY = (int) Math.ceil(d.MAX_Y/d.BOX_SIDE);

		int  xDraw(int x,int sub){return (int) (x*scale)+margin-sub;}
		int  yDraw(int y, int sub){return(int)((d.MAX_Y-y)*scale)+ margin-sub;}
		double xDdraw(double x, int sub){return (double) (x*scale)+margin-sub;}
		double yDdraw(double y, int sub){return (double) ((d.MAX_Y-y)*scale)+margin-sub;}

		 public void paintComponent(Graphics g) {
			super.paintComponent(g);
			 g.setColor(Color.black);
			 // draw grid, first horzontal, then vertical
			 for (int i = 0; i < d.n; i++){
			    drawPoint(i,g);
		     }
             int i;

        //   System.out.println(" XM:"+XM+", YM:"+YM+"restX:"+restX+
         //                    ", restY:"+restY+", d.BOX_SIDE:"+d.BOX_SIDE);
		     g.setColor(Color.green);
		     for (i = 0 ; i <= (numX+1)*d.BOX_SIDE; i += d.BOX_SIDE) {
		       g.drawLine(xDraw(i,0), yDraw(0,0),xDraw(i,0), yDraw( (numY+1)*d.BOX_SIDE,0));
		     }

		      for ( i = 0; i <= (numY+1)*d.BOX_SIDE; i += d.BOX_SIDE) {
			 	 g.drawLine(xDraw(0,0), yDraw(i,0),xDraw((numX+1)*d.BOX_SIDE,0),  yDraw(i,0));
		     }

			  g.setColor(Color.red);
			 // draw cohull
			 int x2 = d.x[d.theCoHull.dtEdges[0]], y2 = d.y[d.theCoHull.dtEdges[0]],x1,y1;
			 for ( i = 1; i < d.theCoHull.numElem; i++){
				 y1 = y2; x1=x2;
				 x2 = d.x[d.theCoHull.dtEdges[i]];
				 y2 = d.y[d.theCoHull.dtEdges[i]];
		         g.drawLine (xDraw(x1,0),yDraw(y1,0), xDraw(x2,0),yDraw(y2,0));
			 }


			  g.drawLine (xDraw(d.x[d.theCoHull.dtEdges[d.theCoHull.numElem-1]],0),
			              yDraw(d.y[d.theCoHull.dtEdges[d.theCoHull.numElem-1]],0),
		                  xDraw(d.x[d.theCoHull.dtEdges[0]],0),yDraw(d.y[d.theCoHull.dtEdges[0]],0));

			 // Draw Delaunay Edges
		     for (int p = 0; p < d.n ; p++) {
				 for ( i = 0; i <= d.delEdges[p].length; i++) {
					 int q =0;
					 if( i== d.delEdges[p].length) q = d.delEdges[p][0];else q =d.delEdges[p][i];
					 // for every DE from 'p'
					 g.drawLine (xDraw(d.x[p],0),yDraw(d.y[p],0), xDraw(d.x[q],0),yDraw(d.y[q],0));
			     }
				 }
			// Draw 10m kote
				int[] naboer = d.beregnMuligeKoteNaboer(d.delEdges[0][0], d.delEdges[0][1]);
				System.out.println("[0][0] = "+d.delEdges[0][0] + " [0][1] = " + d.delEdges[0][1]);
				//for(i = 0; i < naboer.length; i++){
					double[] p1 = d.kotetrekking(d.delEdges[0][0], naboer[0], 40);
					System.out.println(d.delEdges[0][0] + " " + naboer[1]);
					double[] p2 = d.kotetrekking(d.delEdges[0][0], naboer[1], 40);
					//System.out.println("p1_0 "+ p1[0] + "p1_1 "+ p1[1] +"p2_0 "+ p2[0] +"p2_1 "+ p2[1]);
					System.out.println(p2[0] +", "+ p2[1] + " to "+ p1[0] + ", "+ p1[1]);
					Graphics2D g2 = (Graphics2D) g;
					//g2.draw(new Line2D.Double(xDraw(1,0), yDraw(3,0), xDraw(1,0), yDraw(1,0)));
					g2.draw(new Line2D.Double(xDdraw(p1[0],0) , yDdraw(p1[1],0), xDdraw(p2[0],0),yDdraw(p2[1],0)));
					
				//}


		  } // end paintComponent

		  /** called from ge */
			 public void tegnUt() { repaint();}

	}  // end class Graph


   void drawDT(String s){
	 setTitle(s);
	 grafen.tegnUt();
	// Temporary dump
	d.println("drawDT: "+ s);
	int max = Math.min (d.n, 100);
	for (int i = 0;i < max; i++){
		d.println("x["+i+"]="+d.x[i]+ ", y["+i+"]="+d.y[i] + ", z["+i+"]="+d.z[i]);
	  }
   }// end drawDT

}// end class DT
