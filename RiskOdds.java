import java.util.*;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;

public class RiskOdds {
    public static void main(String[] args) {
        
        System.out.println("Hello Cruel World!");
        
        int A = 10;
        int D = 10;
        
        transientStatesMatrix transientStates = new transientStatesMatrix(A,D);
        DoubleMatrix Q = transientStates.matrix();
//        transientStates.print();
        
        absorbingStatesMatrix absorbingStates = new absorbingStatesMatrix(A,D);
        DoubleMatrix R = absorbingStates.matrix();
//        absorbingStates.print();

        DoubleMatrix I = DoubleMatrix.eye(A*D);
        riskMatrix unityMatrix = new riskMatrix(I);
//        unityMatrix.print();

        DoubleMatrix N = new DoubleMatrix(A*D, A*D);
        N = Solve.pinv(I.sub(Q));
        riskMatrix fundamentalMatrix = new riskMatrix(N);
//        fundamentalMatrix.print();
        
        DoubleMatrix F = new DoubleMatrix(A*D, A+D);
        F = N.mmul(R);
        riskMatrix resultMatrix = new riskMatrix(F);
        resultMatrix.print();
    }
    
    public static class transientStatesMatrix extends transitionMatrix {

        // constructor
        public transientStatesMatrix(int attackingArmies, int defendingArmies) {
            // call the parent constructor
            super(attackingArmies, defendingArmies);
            
            // the matrix is armiesA*armiesD x armiesA*armiesD
            int dimensions = armiesA * armiesD;
            
            // instantiate the matrix and populate initially with zeros
            matrix = DoubleMatrix.zeros(dimensions,dimensions);
            
            // populate the matrix with the probabilities
            populateMatrix();
        }
        
        // we need our own populateRow() method for transient states
        protected DoubleMatrix populateRow(int row) {
            DoubleMatrix thisRow = DoubleMatrix.zeros(1,matrix.columns);
            
            // how many armies the attacker/defender has for this row/state
            int aIn = row % armiesA + 1;
            int dIn = row / armiesA + 1;
            
            //            String message = "Row " + row + ": (" + aIn + "," + dIn + ") -- ";
            
            // loop through and populate odds for each column of this row
            // we don't need to go further out than the diagonal, because those values are always zero
            for (int col=0; col<row; col++) {
                // how many armies the attacker and defender have for each output state
                int aOut = col % armiesA + 1;
                int dOut = col / armiesA + 1;
                
                //                message += " (" + aOut + "," + dOut + ")";
                
                // how many dice were rolled this turn
                int aDice = Math.min(aIn,3);
                int dDice = Math.min(dIn,2);
                int aDiff = aIn - aOut;
                int dDiff = dIn - dOut;
                double odds = getOdds(aDice, dDice, aDiff, dDiff);
                
                thisRow.put(col, odds);
            }
            
            //            System.out.println(message);
            
            return thisRow;
        }
    }
    
    public static class absorbingStatesMatrix extends transitionMatrix {
        
        // constructor
        public absorbingStatesMatrix(int attackingArmies, int defendingArmies) {
            // call the parent constructor
            super(attackingArmies, defendingArmies);
            
            // instantiate the matrix and populate initially with zeros
            // the matrix is armiesA*armiesD x armiesA + armiesD
            matrix = DoubleMatrix.zeros(armiesA * armiesD,armiesA + armiesD);
            
            // populate the matrix with the probabilities
            populateMatrix();
        }
        
        // we need our own populateRow() method for absorbing states
        protected DoubleMatrix populateRow(int row) {
            DoubleMatrix thisRow = DoubleMatrix.zeros(1,matrix.columns);
            
            // how many armies the attacker/defender has for this row/state
            int aIn = row % armiesA + 1;
            int dIn = row / armiesA + 1;
            
//            String message = "Row " + row + ": (" + aIn + "," + dIn + ") -- ";
            
            // loop through and populate odds for each column of this row
            for (int col=0; col<thisRow.columns; col++) {
                // how many armies the attacker and defender have for each output state
                int aOut = Math.max(0, col - armiesD + 1); //col % armiesA + 1;
                int dOut = col < armiesD ? col + 1 : 0; //col / armiesA + 1;
                
//                message += " (" + aOut + "," + dOut + ")";
                
                // how many dice were rolled this turn
                int aDice = Math.min(aIn,3);
                int dDice = Math.min(dIn,2);
                int aDiff = aIn - aOut;
                int dDiff = dIn - dOut;
                double odds = getOdds(aDice, dDice, aDiff, dDiff);
                
                thisRow.put(col, odds);
            }
            
//            System.out.println(message);
            
            return thisRow;
        }
    }
    
    // transientMatrix class forms a transition matrix
    protected static class transitionMatrix extends riskMatrix {
        protected int armiesA;
        protected int armiesD;
        
        // constructor
        public transitionMatrix(int attackingArmies, int defendingArmies) {
            super();
            this.armiesA = attackingArmies;
            this.armiesD = defendingArmies;
        }
        
        // override parent print method so we can display the number of armies for the matrix
        public void print() {
            System.out.println(armiesA + " attacking armies, " + armiesD + " defending armies");
            super.print();
        }
        
        // protected method to populate the matrix with probabilities
        protected void populateMatrix() {
            // loop over whole matrix and populate it with the probabilities
            for(int i=0; i<matrix.rows; i++) {
                matrix.putRow(i,populateRow(i));
            }
        }
        
        // populate a given row of the matrix
        protected DoubleMatrix populateRow(int row) {
            DoubleMatrix thisRow = DoubleMatrix.zeros(1,matrix.columns);
            
            
            return thisRow;
        }
    
        // retrieve the odds of the given losses with the given number of dice
        protected double getOdds(int aDice, int dDice, int aDiff, int dDiff) {
            int totalLoss = Math.min(aDice,dDice);
            double odds = 0.0d;
            
            if (aDiff >= 0 && aDiff <= totalLoss && dDiff >= 0 && dDiff <= totalLoss && aDiff + dDiff == totalLoss) {
                
                // 3 v 2
                if (aDice == 3 && dDice == 2) {
                    // a lose 2
                    if (aDiff == 2) {
                        odds = 2275d / 7776d;
                    }
                    // split
                    if (aDiff == 1) {
                        odds = 2611d / 7776d;
                    }
                    // d lose 2
                    if (dDiff == 2) {
                        odds = 2890d / 7776d;
                    }
                }
                
                // 3 v 1
                if (aDice == 3 && dDice == 1) {
                    // a loses 1
                    if (aDiff == 1) {
                        odds = 441d / 1296d;
                    }
                    // d loses 1
                    if (dDiff == 1) {
                        odds = 855d / 1296d;
                    }
                }
                
                // 2 v 2
                if (aDice == 2 && dDice == 2) {
                    // a loses 2
                    if (aDiff == 2) {
                        odds = 581d / 1296d;
                    }
                    // split
                    if (aDiff == 1) {
                        odds = 420d / 1296d;
                    }
                    // d loses 2
                    if (dDiff == 2) {
                        odds = 295d / 1296d;
                    }
                }
                
                // 2 v 1
                if (aDice == 2 && dDice == 1) {
                    // a loses 1
                    if (aDiff == 1) {
                        odds = 91d / 216d;
                    }
                    // d loses 1
                    if (dDiff == 1) {
                        odds = 125d / 216d;
                    }
                }
                
                // 1 v 2
                if (aDice == 1 && dDice == 2) {
                    // a loses 1
                    if (aDiff == 1) {
                        odds = 161d / 216d;
                    }
                    // d loses 1
                    if (dDiff == 1) {
                        odds = 55d / 216d;
                    }
                }
                
                // 1 v 1
                if (aDice == 1 && dDice == 1) {
                    // a loses 1
                    if (aDiff == 1) {
                        odds = 21d / 36d;
                    }
                    // d loses 1
                    if (dDiff == 1) {
                        odds = 15d / 36d;
                    }
                }
            }
            
            return odds;
        }
    }
    
    public static class riskMatrix {
        protected DoubleMatrix matrix;
        
        // constructor
        public riskMatrix() {
            
        }
        // and a version with the matrix supplied
        public riskMatrix(DoubleMatrix m) {
            matrix = m;
        }

        // method that returns the actual matrix
        public DoubleMatrix matrix() {
            return matrix;
        }
        
        // public method to print the matrix
        public void print() {
            for (int row=0; row<matrix.rows; row++) {
                String thisRow = "";
                for (int col=0; col<matrix.columns; col++) {
                    double el = matrix.get(row,col) + 0.0d; // (we add 0.0 to get rid of negative signed zeros)
                    String elStr = el < 0 ? "" : " "; // add a space for positive numbers to line up with any negative signs
                    elStr += String.format("%.2f", el);
                    thisRow += " " + elStr;
                }
                System.out.println(thisRow);
            }
        }
    }
}