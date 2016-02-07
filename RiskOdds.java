import java.util.*;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;

public class RiskOdds {
    public static void main(String[] args) {
        
        System.out.println("CURRENTLY BROKEN, SEE CODE AT LINE 483");
        
        int A = -1; // number of attacking armies
        int D = -1; // number of defending armies
        int places = -1; // number of decimal places to display in table
        
        // for single battles, we use the actual dice odds depending on the number of armies in each roll
        // but if we're calculating the odds of a battle in the middle of a path,
        // we need to find not the odds of winning the battle, but the odds of winning the battle
        // with a bunch of armies left over (however many extra we want to save for the rest of the path)
        // so the odds will assume we're always using 3 dice
        // so the following boolean should be false if we want to find the odds of a intermediary battle
        boolean isTerminalBattle = false;
        
        // get inputs for A, D, and places
        Scanner in = new Scanner(System.in).useDelimiter("\n");
        do {
            System.out.print("Enter the number of Attacking armies (integer >= 1): ");
            if (in.hasNextInt()) { A = in.nextInt(); }
            else { in.next(); }
        } while (A < 1);
        do {
            System.out.print("Enter the number of Defending armies (integer >= 1): ");
            if (in.hasNextInt()) { D = in.nextInt(); }
            else { in.next(); }
        } while (D < 1);
        do {
            System.out.print("Enter the number of decimal places to display (0 <= integer <= 10): ");
            if (in.hasNextInt()) { places = in.nextInt(); }
            else { in.next(); }
        } while (places < 0 || places > 10);

        // create transient states matrix 'Q'
        transientStatesMatrix transientStates = new transientStatesMatrix(A,D, isTerminalBattle);
        DoubleMatrix Q = transientStates.matrix();
//        transientStates.print();
        
        // create absorbing states matrix 'R'
        absorbingStatesMatrix absorbingStates = new absorbingStatesMatrix(A,D, isTerminalBattle);
        DoubleMatrix R = absorbingStates.matrix();
//        absorbingStates.print();

        // create unity matrix 'I'
        DoubleMatrix I = DoubleMatrix.eye(A*D);
        riskMatrix unityMatrix = new riskMatrix(I);
//        unityMatrix.print();

        // create matrix 'N' == (I-Q)^-1
        DoubleMatrix N = new DoubleMatrix(A*D, A*D);
        N = Solve.pinv(I.sub(Q));
        riskMatrix fundamentalMatrix = new riskMatrix(N);
//        fundamentalMatrix.print();
        
        // create final probability matrix 'F' == N * R
        DoubleMatrix F = new DoubleMatrix(A*D, A+D);
        F = N.mmul(R);
        riskMatrix resultMatrix = new riskMatrix(F);
//        resultMatrix.print();
        
        // create table of probabilities for attacker winning each matchup up to A vs D
        double threshold = 0.78d;
        probabilitiesTable probTable = new probabilitiesTable(A,D,F);
        probTable.print(places);
        probTable.printArmiesNeeded(threshold);
        
        probTable.temp(threshold);
    }
 
    // creates a matrix of the probability of the attacker winning the battle,
    // taking the form of a table in which each column represents a different number of attacking armies
    // in ascending order, from 1 to 'attackingArmies'
    // and each row represents a different number of defending armies from 1 to 'defendingArmies',
    // and each cell represents the odds of the attacker winning that matchup
    public static class probabilitiesTable {
        protected int armiesA;
        protected int armiesD;
        protected DoubleMatrix matrix;
       
        // constructor
        public probabilitiesTable(int attackingArmies, int defendingArmies, DoubleMatrix probMatrix) {
            this.armiesA = attackingArmies;
            this.armiesD = defendingArmies;
            
            this.matrix = constructTable(probMatrix);
        }
        
        // method that returns the table as the actual matrix
        public DoubleMatrix matrix() {
            return matrix;
        }
        
        // display table with labeled axes to precision 'places'
        public void print(int places, double threshold) {
            // set formatting for precision
            if (places < 0) {
                places = 0;
            }
            String placesString = "%." + places + "f";
            
            // print column labels at the top (attacking armies)
            String topLabel = spaces(digits(armiesD)) + "  ";
            String line = topLabel;
            for (int col=0; col<matrix.columns; col++) { // loop through columns to add number
                topLabel += " " + (col + 1) + spaces(places + 2 - digits(col + 1)) + " "; // add column number (== attacking armies) + spacing
            }
            for (int i=line.length(); i<topLabel.length(); i++) { // add appropriate number of underscores to 'line'
                line += "_";
            }
            System.out.println(topLabel);
            System.out.println(line);
            
            // print values (along with row labels (defending armies))
            for (int row=0; row<matrix.rows; row++) {
                String thisRow = "" + (row + 1) + spaces(digits(row + 1, armiesD)) + " |"; // the row label is the number of defending armies (plus some spaces for alignment purposes)
                for (int col=0; col<matrix.columns; col++) {
                    double el = matrix.get(row,col) + 0.0d; // (we add 0.0 to get rid of negative signed zeros)
                    if (el < threshold) {
                        thisRow += " " + spaces(places + 2) + " ";
                    } else {
                        thisRow += " " + String.format(placesString, el) + " "; // display to 2 decimal places of precision
                    }
                }
                System.out.println(thisRow);
            }
        }
        public void print(int places) {
            print(places, -1000d);
        }
        public void print() {
            print(2); // print with 2 decimal places by default
        }
        
        // return the difference between the number of digits of 'maxNum' and 'thisNum'
        private int digits(int thisNum, int maxNum) {
            // find the difference in digits between maxNum and thisNum
            int numDigits;
            if (thisNum < 0) { // if thisNum is less than 0, treat it as 0 digits long
                numDigits = (int) Math.log10(maxNum) + 1;
            } else {
                numDigits = (int) Math.log10(maxNum) - (int) Math.log10(thisNum);
            }
            return numDigits;
        }
        // overloaded single parameter version, pass -1 for thisNum to ignore it
        private int digits(int num) {
            return digits(-1, num);
        }
        
        // return <digits> spaces as a string, for display alignment purposes
        private String spaces(int digits) {
            String spaces = "";
            for (int i=0; i<digits; i++) {
                spaces += " ";
            }
            return spaces;
        }
        
        // take the probabilities matrix and return the table (also a matrix)
        protected DoubleMatrix constructTable(DoubleMatrix probMatrix) {
            DoubleMatrix table = new DoubleMatrix(armiesD, armiesA);
            
            for (int row=0; row<table.rows; row++) {
                for (int col=0; col<table.columns; col++) {
                    double prob = 0.0d;
                    // loop through probMatrix for the appropriate row and get the sum
                    // of all probabilities in absorbing states where the attacker wins
                    // (i.e. the last 'armiesA' columns of that row)
                    // this is the probability that the attacker wins for this matchup
                    int i = row * armiesA + col; // i is the row to look up in probMatrix
                    for (int j=armiesD; j<armiesA+armiesD; j++) { // j are the columns we want to add up in probMatrix
                        prob += probMatrix.get(i,j);
                    }
                    // set the probability for this matchup
                    table.put(row,col,prob);
                }
            }
            
            return table;
        }
        
        // print a number of defending armies + number of attacking armies it will take to beat them at the given odds
        public void printArmiesNeeded(int defenders, double winningOdds) {
            for (int col=0; col<matrix.columns; col++) {
                if (matrix.get(defenders-1,col) >= winningOdds) {
                    System.out.println(defenders + " " + (col + 1));
                    break;
                }
            }
        }
        public void printArmiesNeeded(double winningOdds) {
            for (int row=0; row<matrix.rows; row++) {
                printArmiesNeeded(row+1, winningOdds);
            }
        }
        
        public void temp(double winningOdds) {
            for (int row=0; row<matrix.rows; row++) {
                int defenders = row + 1;
                for (int col=0; col<matrix.columns; col++) {
                    if (matrix.get(defenders-1,col) >= winningOdds) {
                        int attackers = col + 1;
                        double padding = attackers - ((double) defenders * (7161d / 8391d));
                        System.out.println(defenders + " " + padding);
                        break;
                    }
                }
            }
        }
    }
    
    public static class transientStatesMatrix extends transitionMatrix {

        // constructor
        public transientStatesMatrix(int attackingArmies, int defendingArmies, boolean isTerminalBattle) {
            // call the parent constructor
            super(attackingArmies, defendingArmies, isTerminalBattle);
            
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
        public absorbingStatesMatrix(int attackingArmies, int defendingArmies, boolean isTerminalBattle) {
            // call the parent constructor
            super(attackingArmies, defendingArmies, isTerminalBattle);
            
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
    
    // transitionMatrix class forms a transition matrix
    protected static class transitionMatrix extends riskMatrix {
        protected int armiesA;
        protected int armiesD;
        protected boolean isTerminal;
        
        // constructor
        public transitionMatrix(int attackingArmies, int defendingArmies, boolean isTerminalBattle) {
            super();
            this.armiesA = attackingArmies;
            this.armiesD = defendingArmies;
            this.isTerminal = isTerminalBattle;
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
                
                // isTerminal is a boolean which tells us whether we're calculating for a regular battle
                // or an intermediary battle. If it's true, we want to use the regular odds.
                // if it's false, then we want to calculate for an intermediary battle,
                // where we want to know the odds of winning the battle as though we always got to roll
                // 3 dice. That's because in the game, we'll have extra armies on the country that we're
                // saving for future battles along the path
                if (isTerminal == true) {
                
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
                } else {
                    // if we're here, isTerminal is false, so we want to pretend the attacker
                    // always gets to roll three dice, no matter how few armies he has
                    
                    // 3 v 2 and 2 v 2
                    if ((aDice == 3 || aDice == 2) && dDice == 2) {
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
                    
                    // 3 v 1 and 2 v 1 and 1 v 1
                    if (dDice == 1) {
                        // a loses 1
                        if (aDiff == 1) {
                            odds = 441d / 1296d;
                        }
                        // d loses 1
                        if (dDiff == 1) {
                            odds = 855d / 1296d;
                        }
                    }
                    
                    // 1 v 2
                    
                    // THIS BLOCK IS CURRENTLY BROKEN;
                    // NEED TO FIGURE OUT THE ACTUAL STATES THE BATTLE CAN GO INTO
                    // WHEN ROLLING 3 DICE WITH ONLY ONE ARMY
                    // SOME OF THOSE STATES WILL VIOLATE THE LONG CONDITIONAL
                    // WE SET UP ABOVE, SO WE NEED TO MOVE OUT OF IT?
                    
                    if (aDice == 1 && dDice == 2) {
                        // a loses 1
                        if (aDiff == 1) {
                            // (these are actually the 3v2 odds of a split plus the odds of attacker losing 2)
                            odds = 4886d / 7776d;
                        }
                        // d loses 1
                        if (dDiff == 1) {
                            // (these are actually the 3v2 odds of defender losing 2)
                            odds = 2890d / 7776d;
                        }
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
                    elStr += String.format("%.2f", el); // display to 2 decimal places of precision
                    thisRow += " " + elStr;
                }
                System.out.println(thisRow);
            }
        }
    }
}