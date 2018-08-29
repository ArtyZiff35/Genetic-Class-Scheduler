/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package oma_assignment;

import java.awt.AWTEventMulticaster;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import java.util.Arrays;

/**
 * @author Arturo
 */

public class OMA_Assignment_1 {
    
    private static int NumTimeSlot = -1;
    private static int NumEsami = -1;
    private static int NumStudenti = 0;
    private static boolean[][] StudentExamMatrix;
    private static int[][] CommonStudentMatrix;
    private static int Iteration = 0;
    private static Exam[] ExamArray;
    private static int[] Solution;
    private static int[] TempSolution;
    private static float Penalty = -1;
    private static float TempPenalty = -1;
    private static int[] RowMoveRecord = new int[3];
    private static boolean GreenFLAG=true;
    
    private static String InstanceName="";
    private static String TimeName="";

    public static void main(String[] args) throws FileNotFoundException, IOException {
       
        if(args.length!=3)
        {
            System.err.println("Wrong number of arguments!");
            return;
        }
        else
        {
            InstanceName=args[0];
            TimeName=args[2];
        }
        
        //PHASE 0: instance elaboraion  //////////////////////////////////////////////////////////////////////////////////
        ReadFromFile();                         //This method simply reads the 3 input files and builds the StudentExamMatrix, other than calculating NumTimeSlot, NumEsami and NumStudenti
        CalculateCommonStudentMatrix();         //This method calculates the number of common students between each couple of exams (CommonStudentMatrix)
        Calculate_Exam_Array();                 //Here we create the array of Exams (ExamArray), where for each exam we store its total number of conflicts and the number of total common students with other exams
        //PHASE 1: generate initial feasable solution ////////////////////////////////////////////////////////////////////
        
        int FileCounter=1;
        
        
        Solution = new int[NumEsami];           //We initialize the final solution
        for(int i=0; i<Solution.length; i++)
            Solution[i]=-1;
        
        //PHASE 2: metaheuristics   ///////////////////////////////////////////////////////////////////////////////////
        //If IterationForSolution is -1, that means no initial solution was found!
        
        int[] FinalSolution = new int[NumEsami];
        float FinalPenalty = 1000f;
        
        
        /* The following commented lines, including the if, have to be used if we want to execute
         * one of the other MetaHeuristic approaches
         */
        
        //GenerateInitialSolution(300, true);    //We generate a feasable solution using tot max iterations. If a solution here was found, it WILL CERTAINLY be feasible, no need to check afterwards
        //if(Solution != null)
        //{
        /*System.out.println("Initial Solution found with penalty " + Penalty + " at iteration " + Iteration);*/
            
        new java.util.Timer().schedule( 
                    new java.util.TimerTask() {
                        @Override
                        public void run() {
                            // Code to be executed at a specific time...
                            System.out.println("STOP");
                            GreenFLAG = false;                            
                        }
                    }, 
                    (Integer.parseInt(TimeName)*1000) //Seconds * 1000
            );
            FinalSolution = GeneticAlgorithm();
            FinalPenalty = CheckPenaltyFinalSolution(FinalSolution);
            WriteFile(FinalSolution,FinalPenalty,FileCounter++);
            GreenFLAG = true;
            
        //}
        
       
    }
    
    private static int[] GeneticAlgorithm()
    {
        //We generate a random population
        Random r = new Random();
        int[] BestSoFar = new int[NumEsami];
        int PopulationSize = 300;
        int[][] Population = new int[PopulationSize][NumEsami];
        float[] PopulationPenalty = new float[PopulationSize];
        for(int i=0; i<PopulationSize; i++)
        {
            Population[i] = GenerateInitialSolution(300, false);  
            PopulationPenalty[i] = CheckPenaltyFinalSolution(Population[i]);
        }
        System.out.println("Starting GENETIC with a Population of " + PopulationSize );
        int GeneticCounter = 0;
        int LastBestIteration = 0;
        int REDZONE = 10000;
        boolean Scramble = false;
        float CurrentGlobalMIN = 9999999;
        float SogliaMutazione = 0.30f;
        float SogliaSwap = 0.8f;
        float SogliaShuffle = 0.03f;
        float SogliaCrossOver = 0.9f;
        //We start the loop
        while(GreenFLAG == true)
        {
            GeneticCounter++;
            //FIRST, we select some solutions, and two of them become parents (the best two)
            int[] Parent1 = new int[NumEsami];
            int[] Parent2 = new int[NumEsami];
            int[] Child = new int[NumEsami];
            int PotentialSize = 4;
            int[] PotentialParentsIndexes = new int[PotentialSize];
            for(int i=0; i<PotentialSize; i++)
                PotentialParentsIndexes[i] = -1;    //We need to initialize this array to see if it will be full
            //Here we find the POTENTIAL parents
            int counter = 0;
            while(counter < PotentialSize)
            {
                boolean found = true;
                int TmpIndex = -1;
                while(found == true)
                {
                    //We generate a random population index
                    TmpIndex = ThreadLocalRandom.current().nextInt(0,PopulationSize);
                    //We see if it was already chosen in this iteration
                    for(int i=0; i<PotentialSize; i++)
                    {
                        if(PotentialParentsIndexes[i] != TmpIndex)
                            found = false;
                    }
                }
                //We increment the counter and save the index we found as a potential parent
                PotentialParentsIndexes[counter] = TmpIndex;
                counter++;      
            }
            
            //BUBBLE SORT to sort potential parents
            for(int i = 0; i < PotentialParentsIndexes.length; i++)
            {
                boolean flag = false;
                for(int j = 0; j < PotentialParentsIndexes.length-1; j++) {

                    if(PopulationPenalty[PotentialParentsIndexes[j]]>PopulationPenalty[PotentialParentsIndexes[j+1]]) {
                        int k = PotentialParentsIndexes[j];
                        PotentialParentsIndexes[j] = PotentialParentsIndexes[j+1];
                        PotentialParentsIndexes[j+1] = k;
                        flag=true; 
                    }


                }

                if(!flag) break; 
            }
            Parent1 = Population[PotentialParentsIndexes[0]];
            Parent2 = Population[PotentialParentsIndexes[1]];
            
            
            
            //Now we do a UNIFORM CROSS OVER between them and generate a CHILD
            float CrossChance = r.nextFloat();
            if(CrossChance < SogliaCrossOver)
                Child = CrossOver(Parent1,Parent2);
            else
                Child = Parent1;
            
              
            
            //++++++++EVENTUAL CHILD MUTATIONS SHOULD BE DONE HERE+++++++++++
            float Chance = r.nextFloat();
            if(Chance < SogliaMutazione)
            {
                float NewChance = r.nextFloat();
                if(NewChance<0.4f)
                    Child = GenerateRandomNeighbor(Child,1);    
                else if(NewChance < 0.6f)
                    Child = SwapTwoDays(Child);            
                else if(NewChance < 0.65f)
                    Child = ShuffleArray(Child);
                else
                {
                    //First mutation is moving one or more exam 
                    Child = GenerateRandomNeighbor(Child,1);
                    //Second mutation is swapping two full time slots
                    float SwapChance = r.nextFloat();
                    if(SwapChance <= SogliaSwap)
                        Child = SwapTwoDays(Child);            
                    //Third mutation is using Gabri's algorithm
                    Chance = r.nextFloat();
                    if(Chance < SogliaShuffle)
                        Child = ShuffleArray(Child);
                }
            } 
                 
            
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                //Stagnation is detected: we try to make bigger jumps by changing chance of mutation
            if(Scramble==true)
            {
                System.out.println("Stagnation was detected. Elevating Chance limits...");
                Scramble=false;
                SogliaMutazione += 0.1f;
                SogliaSwap += 0.05f;
                SogliaShuffle += 0.02f;
                LastBestIteration = GeneticCounter;   
                REDZONE += 20000;
            }
            
            //++++++++++++++++++++++++++++++++++++++++++++++++  
                
            //Now we check if that Child already exists in the population, in that case we DO NOT add him to the population
            if(CheckChildExistance(Child,Population,PopulationSize)==false)
            {
                //We put the child in place of the last (less optimized) population member who was a potential parent
                int index;
                index = PotentialParentsIndexes[PotentialParentsIndexes.length-1];
                Population[index] = Arrays.copyOf(Child, NumEsami);
                PopulationPenalty[index] = CheckPenaltyFinalSolution(Child);
                
                ////// In case we find a new minimum/////////////////////
                if(PopulationPenalty[index] < CurrentGlobalMIN)
                {
                    CurrentGlobalMIN = PopulationPenalty[index];
                    BestSoFar = Child;
                    
                    //Iterations stuff...
                    LastBestIteration = GeneticCounter;
                                        
                    System.out.println("Current best : " + CurrentGlobalMIN + "  at iteration " + LastBestIteration);
                }
            }
            
            //If we are stuck, modify something...
            if((GeneticCounter-LastBestIteration) >= REDZONE)
                Scramble = true;
            
        }
        
        return BestSoFar;
 
    }
    
    private static int[] CrossOver(int[] Parent1, int[] Parent2)
    {
        int[] Child = new int[NumEsami];
        int[] OtherParent = new int[NumEsami];
        //We swap semi-casually each exam's time slot with the one of the other parent
        
        //We start randomically from one of the two
        
        Random r = new Random();
        float StartingChance = r.nextFloat();
        if(StartingChance <= 0.5f)  
        {
            Child = Arrays.copyOf(Parent1, NumEsami);;
            OtherParent = Arrays.copyOf(Parent2, NumEsami);;
        }
        else
        {
            Child = Arrays.copyOf(Parent2, NumEsami);;
            OtherParent = Arrays.copyOf(Parent1, NumEsami);;
        }
        
            
        
        //Now we actually iterate through the Child's exams and probablistivly swap their timeslot
        for(int i=0; i<NumEsami; i++)
        {
            float SwappingChance = r.nextFloat();
            if(SwappingChance <= 0.5f)
            {
                int BackIndex = Child[i];   //In case we have to switch back
                Child[i] = OtherParent[i];
                ArrayList<Integer>[] ConvertedChild = ReconvertSolution(Arrays.copyOf(Child, NumEsami));
                if(CheckConflictInTempSolution(ConvertedChild[Child[i]],i)==false)
                {
                    //The current solution was NOT feasible, so we turn back this move
                    Child[i] = BackIndex;
                }
            }
            
        }
        return Child;
    }
    
    private static int[] GenerateRandomNeighbor(int[] Solution, int num)
    {
        // num value is the number of exams to be taken from the solution and reassigned to a timeslot
        boolean go;
        int[] Temp1;            //This temporary solution contains the modification applied trying to build the neighbor
        int[] Temp2 = Arrays.copyOf(Solution, NumEsami); //This temporary solution instead, contains the corrt partial neighbor solution, at the end the competerd neighbour
                
        for(int n = 0; n < num; n++) 
        {
            go = true;
            while(go)
            {
                // nextInt is normally exclusive of the top value,
                // so add 1 to make it inclusive
                int random1 = ThreadLocalRandom.current().nextInt(0, NumEsami);
                
                Temp1 = Arrays.copyOf(Temp2, NumEsami);
                ArrayList<Integer>[] ConvertedArray = ReconvertSolution(Temp1);

                Exam e = ExamArray[random1];
                e.clearFreeSlots();
                int counter = 0;
                
                //This loop counts and store in the Exam e object the available timeslot
                for(int index = 0; index < ConvertedArray.length; index++)  //ConvertedArray.lenght is the number of TimeSlots
                {
                    if(index != Temp1[random1]) //Make sure we aren't working on the same slot
                    {
                        //IF CELL IS EMPTY
                        if(ConvertedArray[index].isEmpty())
                        {
                            counter++;
                            e.addFreeSlot(index); //the exam can be put inside the TimeSlot index
                        }
                        //IF CELL IS ALREADY FULL BUT NO CONFLICTS 
                        else if(CheckConflictInTempSolution(ConvertedArray[index],e.getExamID())==true)
                        {
                            counter++;
                            e.addFreeSlot(index);
                        }
                    }
                }
                //We save the number of possible free slots into the exam
                e.setFreeSlots(counter);


                //IF NO TIME SLOTS AVAILABLE, RESTART
                if(e.getFreeSlots()==0)
                {
                    //System.out.println("No free slot.");
                }
                //OTHERWISE WE SWAP
                else
                {   
                    Random r = new Random();
                    float chance = r.nextFloat();
                    //We swap for a better solution 50% of the time
                    if(chance <= 0.4f)
                    {
                        //I find the best slot in which the exam can be put, instead of choosing it randomly
                        int[] Temp3;
                        float CurrentPenalty;
                        float BestPenalty=1000f;
                        for(int i:e.returnSlots())
                        {
                            Temp3 = Arrays.copyOf(Temp1, NumEsami);
                            Temp3[random1] = i;
                            CurrentPenalty = CheckPenaltyFinalSolution(Temp3);
                            if(CurrentPenalty<BestPenalty)
                            {
                                BestPenalty = CurrentPenalty;
                                Temp2 = Temp3;
                                RowMoveRecord[0] = Solution[random1];
                                RowMoveRecord[1] = i;
                                RowMoveRecord[2] = random1;
                            }
                            go=false;
                        }
                    }
                    else
                    {
                        int random2 = ThreadLocalRandom.current().nextInt(0, e.getFreeSlots());
                        Temp1[random1] = e.getSlot(random2);
                        
                        go = false;
                        RowMoveRecord[0] = Solution[random1];
                        RowMoveRecord[1] = random2;
                        RowMoveRecord[2] = random1;
                        Temp2 = Arrays.copyOf(Temp1, NumEsami);
                    }
                    
                    
                }
            }
        }
        return Temp2;
    }
    
    private static int[] SwapTwoDays(int[] Child)
    {
        ArrayList<Integer>[] ConvertedChild = ReconvertSolution(Arrays.copyOf(Child, NumEsami));
        ArrayList<Integer> TmpCopy;
        int Day1Index = ThreadLocalRandom.current().nextInt(0,NumTimeSlot);
        int Day2Index = ThreadLocalRandom.current().nextInt(0,NumTimeSlot);
        while(Day1Index == Day2Index)
            Day2Index = ThreadLocalRandom.current().nextInt(0,NumTimeSlot);
        //We swap the content of 2 random timeslots
        TmpCopy = ConvertedChild[Day1Index];
        ConvertedChild[Day1Index] = ConvertedChild[Day2Index];
        ConvertedChild[Day2Index] = TmpCopy;
        //We reconvert the solution back
        int[] Result = ReturnInitialSolutionAsArray(ConvertedChild);
        return Result;
    }
    
    private static int[] SimulatedAnnealingParameterWarmingUp(int NumIterations, int[] TempArray, int LengthNeighborhood, int x)
    {
        int[] Best = new int[NumEsami];
        int Iteration = 0;
        float alpha = 0.99f;
        int L = 10*LengthNeighborhood;
        
        float CurrentMIN = 99999999999f;
        
        System.out.println("Starting SIMU");
        Solution = Arrays.copyOf(TempArray, NumEsami);
        Penalty = CheckPenaltyFinalSolution(Solution);
        
        //Setting Initial temperature
      
        float TPenalty = 0f;
        int[] TSol = new int[NumEsami];
        for(int i=0;i<LengthNeighborhood;i++)
        {
            TSol = GenerateRandomNeighbor(Solution, x);
            TPenalty += CheckPenaltyFinalSolution(TSol);
        }
        TPenalty = TPenalty / LengthNeighborhood;

        float T = (Math.abs(TPenalty - Penalty))/0.69f;

        int counter = 0;
        int stagnation = 1000;
        
        while(counter < NumIterations)
        {
            counter++;
            
            //If the System is stagnating we try to reset the initial temperature , in order to warm it up
            if(stagnation == 0)
            {
                System.out.println("Stagnation Detected");

                //Just recalculate the temperature, from the current solution
                for(int i=0;i<LengthNeighborhood;i++)
                {
                    TSol = GenerateRandomNeighbor(Solution, x);
                    TPenalty += CheckPenaltyFinalSolution(TSol);
                }
                TPenalty = TPenalty / LengthNeighborhood;

                T = (Math.abs(TPenalty - Penalty))/0.69f;
                ///////////////////
                /*float diff = TPenalty-Penalty;
                System.out.println(TPenalty + " "+ Penalty +" " +diff);
                System.out.println(Math.exp(-Math.abs(TPenalty - Penalty) / T));*/
                ///////////////////
             
                stagnation = 500;
            }
            
            ///TEST////
            if(Penalty < CurrentMIN)
            {
                Best = Arrays.copyOf(Solution, NumEsami);
                CurrentMIN = Penalty;
                System.out.println("Current penalty is " + CurrentMIN);
                stagnation = 500;
            }
            else
                stagnation --;
            //////////////////////
            
            TempSolution = GenerateRandomNeighbor(Solution, x);
            TempPenalty = CheckPenaltyFinalSolution(TempSolution);
            if(TempPenalty <= Penalty)
            {
                Solution = Arrays.copyOf(TempSolution, NumEsami);
                Penalty = TempPenalty;
            }
            else
            {                
                float p = (float) Math.exp(-(TempPenalty - Penalty)/T);
                
                Random rand = new Random();
                float CasualP = rand.nextFloat();
                // Check if we have to swap
                if(CasualP < p)
                {
                    Solution = TempSolution;
                    Penalty = TempPenalty;
                }
           
            }
            
            //Prepare for next cycle
            Iteration ++;
            if(Iteration%L==0)
            {
               T=T*alpha; 
            } 
        }
        System.out.println("Best Penalty found " + CurrentMIN);
        return Best;
    }

    private static int[] TabuSearch(int Start[], int NumIterazioni)
    {
        Solution = Start;
        int NeighborhoodDimen = 500;
        float AbsoluteMin = 99999999;
        int TabuLength = 90;
        int[][] Neighborhood = new int[NeighborhoodDimen][NumEsami];
        float [] Penalties = new float[NeighborhoodDimen];
        int[][] TabuList = new int[TabuLength][RowMoveRecord.length];
        int[] CurrentBestSol = Solution;
        
        //LOOP
        int counter=0;
        while(counter < NumIterazioni)
        {
            counter++;
            float CurrentMinValue = 999999999;
            int CurrentMinIndex = -1;
            int[] CurrentTabuRow = new int[RowMoveRecord.length];
            //Now we fill the neighborhood
            for(int i=0; i<NeighborhoodDimen; i++)
            {
                int[] TmpNeigh = GenerateRandomNeighbor(Solution, 1);
                float TmpPen = CheckPenaltyFinalSolution(TmpNeigh);

                //We start to look for the minimum value of the penalty and its index
                if(TmpPen <= CurrentMinValue)
                {
                    CurrentBestSol = TmpNeigh;
                    CurrentMinValue = TmpPen;
                    CurrentTabuRow = RowMoveRecord;
                }

            }
            //If we choose the solution..
            if(CheckConflictInTabuList(CurrentTabuRow,TabuList,TabuLength,CurrentMinValue,AbsoluteMin))
            {
                TabuList = ShiftRowsMatrix(TabuList,TabuLength);
                TabuList[0] = CurrentTabuRow; //We always insert into the first position
                Solution = CurrentBestSol;
                Penalty = CurrentMinValue;
                if(Penalty<=AbsoluteMin) 
                {
                    AbsoluteMin = Penalty;
                    System.out.println("Current penalty in TabuSearch is " + Penalty);
                }
                    
            }
            else
            {
                //Also in this case we save into Solution in order to try to move the solution search forward and not remain stuck into the same neighborhood
                Solution = CurrentBestSol;
            }
        }
        return Solution;
    }
    
    public static void VNS()
    {
        float BEST_ABSOLUTE_PENALTY = 999999f;
        int MAX_NEIGH = 100;
        int NeighborhoodDimension = 20;
        int[][][] Neighborhood = new int[MAX_NEIGH][NeighborhoodDimension][NumEsami];
        float[] BestNeighborhoodPenalty = new float[MAX_NEIGH];
        
        //We generate the initial solution either using the genetic algorithm best solution or a random initial solution  
        //int[] InitialSol = GenerateInitialSolution(300, false);
        int[] InitialSol = GeneticAlgorithm();
        
        //We fill the neighborhoods
        //For each neighborhood
        for(int i=0; i<MAX_NEIGH; i++)
        {
            //For each member, we fill with a random neighbor
            for(int j=0; j<NeighborhoodDimension; j++)
            {
                Neighborhood[i][j] = GenerateRandomNeighbor(InitialSol,ThreadLocalRandom.current().nextInt(0,NeighborhoodDimension/2));
            }
            //WE set the initial penalty
            BestNeighborhoodPenalty[i] = CheckPenaltyFinalSolution(InitialSol);
            
        }
        
        System.out.println("Initial Population created");
        
        //Main Loop
        boolean Stopped = true;
        int CurrentNeighborhoodIndex = 0;
        int[] LocalSearchResult = new int[NumEsami];
        while(true)
        {
            Stopped = true;
            while(Stopped)
            {
                int[] RandomSolution = Arrays.copyOf(Neighborhood[CurrentNeighborhoodIndex][ThreadLocalRandom.current().nextInt(0,NeighborhoodDimension)], NumEsami);
                
                //LOCAL SEARCH che ritorna una LocalSearchResult
                //LocalSearchResult = SteepestDescent(RandomSolution,500);
                LocalSearchResult = FirstImprovement(RandomSolution,200);
                //LocalSearchResult = RandomSolution;    Reduced VNS
                
                float LocalSearchPenalty = CheckPenaltyFinalSolution(LocalSearchResult);
                if(LocalSearchPenalty < BestNeighborhoodPenalty[CurrentNeighborhoodIndex])
                {
                    //Se è migliore, rimaniamo in questo neighborhood, ma lo sostituiamo con quello di LocalSearchResult
                    BestNeighborhoodPenalty[CurrentNeighborhoodIndex] = LocalSearchPenalty;
                    //We create new neighborhood
                    for(int j=0; j<NeighborhoodDimension; j++)
                    {
                        Neighborhood[CurrentNeighborhoodIndex][j] = GenerateRandomNeighbor(LocalSearchResult,ThreadLocalRandom.current().nextInt(0,NeighborhoodDimension/2));
                    }
                    
                    if(LocalSearchPenalty < BEST_ABSOLUTE_PENALTY)
                    {
                        BEST_ABSOLUTE_PENALTY = LocalSearchPenalty;
                        System.out.println("New absolute best found with value " + BEST_ABSOLUTE_PENALTY);
                    }
                    
                }
                else
                {
                    //Se LocalSearchResult NON è migliore..
                    Stopped=false;
                }
                
            }
            CurrentNeighborhoodIndex++;
            if(CurrentNeighborhoodIndex>=MAX_NEIGH)
                CurrentNeighborhoodIndex = 0;
        }
           
    }
    
    private static int[] FirstImprovement(int[] Sol, int Dimension)
    {
        boolean go=true;
        int[] BestSolution = Arrays.copyOf(Sol,NumEsami);
        while(go)
        {
            go=false;
            int[][] Neighborhood = new int[Dimension][NumEsami];
            for(int i=0; i<Dimension; i++)
            {
                Neighborhood[i] = GenerateRandomNeighbor(BestSolution,1);
            }
            
            float BestPenalty = CheckPenaltyFinalSolution(BestSolution);
            float CurrentPenalty=0f;
        
            for(int j=0; j<Dimension; j++)
            {
                CurrentPenalty = CheckPenaltyFinalSolution(Neighborhood[j]);
                if(CurrentPenalty<BestPenalty)
                {
                    BestPenalty = CurrentPenalty;
                    BestSolution = Arrays.copyOf(Neighborhood[j], NumEsami);
                    go=true;
                    j=Dimension;
                }
            }
        }
        
        
        return BestSolution;
    }
    
    private static int[] SteepestDescent(int[] Sol, int Dimension)
    {
        boolean go=true;
        int[] BestSolution = Arrays.copyOf(Sol,NumEsami);
        while(go)
        {
            go=false;
            int[][] Neighborhood = new int[Dimension][NumEsami];
            for(int i=0; i<Dimension; i++)
            {
                Neighborhood[i] = GenerateRandomNeighbor(BestSolution,1);
            }
            
            float BestPenalty = CheckPenaltyFinalSolution(BestSolution);
            float CurrentPenalty=0f;
        
            for(int j=0; j<Dimension; j++)
            {
                CurrentPenalty = CheckPenaltyFinalSolution(Neighborhood[j]);
                if(CurrentPenalty<BestPenalty)
                {
                    BestPenalty = CurrentPenalty;
                    BestSolution = Arrays.copyOf(Neighborhood[j], NumEsami);
                    go=true;
                }
            }
   
        }
        return BestSolution;
    }
    
    public static void WriteFile(int[] temp, float tpen, int counter) throws IOException 
    {
      int n=0,m=0;
      File fout = new File(InstanceName+"_OMAAL_group08.sol");
      FileOutputStream fos = new FileOutputStream(fout);

      BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
      
      for(int i=0; i<NumEsami;i++)
      {
          n=i+1;
          m=temp[i]+1;
          
          if(n<10)
          {
              bw.write("000"+n+" "+m);
              bw.newLine();
          }
          else if(n<100)
          {
              bw.write("00"+n+" "+m);
              bw.newLine();
          }
          else if(n<1000)
          {
              bw.write("0"+n+" "+m);
              bw.newLine();
          }
          else
          {
              bw.write(n+" "+m);
              bw.newLine();
          }
      }

      bw.close();
    }
    
    private static int[] MultiplePointsCrossOver(int[] Parent1, int[] Parent2 , int NumTagli, int NumTentativi)
    {
        int[] Child = new int[NumEsami], TmpParent1 = Arrays.copyOf(Parent1, NumEsami), TmpParent2 = Arrays.copyOf(Parent2, NumEsami);
        int[] Cuts = new int[NumTagli];
        int Counter = NumTentativi;
        boolean Feasibility = false;

        do
        {
            Counter--;
            //We fill with random cuts the array
            for(int i=0; i<NumTagli; i++)
                Cuts[i] = ThreadLocalRandom.current().nextInt(0,NumTimeSlot);
            //We sort the cuts and transform it into an ArrayList
            Arrays.sort(Cuts);
            ArrayList<Integer> CutsList = new ArrayList<Integer>();
            for(int i=0; i<Cuts.length; i++)
                CutsList.add(Cuts[i]);
            //We have a chance of swapping parents
            Random r = new Random();
            float Chance = r.nextFloat();
            if(Chance < 0.5f)
            {
                int[] x = TmpParent1;
                TmpParent1 = TmpParent2;
                TmpParent2 = x;
            }
            //No we actually make the swaps
            boolean FromFirst=true;
            for(int i=0; i<NumEsami; i++)
            {
                //We do the swap first
                if(FromFirst)
                {
                    Child[i] = TmpParent1[i];
                }
                else
                {
                    Child[i] = TmpParent2[i];
                }
                //If we meet a cutting point...
                if(CutsList.contains(i))
                {
                    FromFirst = !FromFirst;
                }
            }
            
            //Finally we check feasibility
            if(Constraint_Check.CheckEverything(Child, CommonStudentMatrix)==true)
            {
                Feasibility = true;
                Counter=-1; //We break the loop
            }

        }while(Counter>0);
        
        if(Feasibility == true)
            return Child;
        else
            return TmpParent1;
        
    }
    
    private static int[] SinglePointCrossOver(int[] Parent1, int[] Parent2 , int NumTentativi)
    {
        int CuttingIndex = ThreadLocalRandom.current().nextInt(1,NumTimeSlot-1); // We exclude the first and last cells
        int[] Child = new int[NumEsami], TmpParent1 = Arrays.copyOf(Parent1, NumEsami), TmpParent2 = Arrays.copyOf(Parent2, NumEsami);
        int Counter = NumTentativi;
        boolean Feasibility = false;
        
        do
        {
            Counter--;
            //First, with 50% side by side chance we copy part of one parent and part of the other
            Random r = new Random();
            float Chance = r.nextFloat();
            if(Chance < 0.5f)
            {
                for(int i=0; i<CuttingIndex; i++)
                    Child[i] = TmpParent1[i];
                for(int i=CuttingIndex; i<NumEsami; i++)
                    Child[i] = TmpParent2[i];
            }
            else
            {
                for(int i=0; i<CuttingIndex; i++)
                    Child[i] = TmpParent2[i];
                for(int i=CuttingIndex; i<NumEsami; i++)
                    Child[i] = TmpParent1[i];
            }
            //Now we simply check the feasibility
            if(Constraint_Check.CheckEverything(Child, CommonStudentMatrix)==true)
            {
                Feasibility = true;
                Counter=-1; //We break the loop
            }
        }while(Counter>0);
        
        //Now we return the Child if it is Feasible, or a Parent Otherwise
        if(Feasibility==true)
            return Child;
        else
            return TmpParent1;
    }
    
    private static int[] ShuffleArray(int[] array)
    {    
        ArrayList<Integer>[] ConvertedTempSolution = ReconvertSolution(Arrays.copyOf(array, NumEsami));
        ArrayList<Integer> temp;
        int index;
        Random random = new Random();
        for (int i = ConvertedTempSolution.length - 1; i >= 0; i--)
        {
            index = random.nextInt(i + 1);
            temp = ConvertedTempSolution[index];
            ConvertedTempSolution[index] = ConvertedTempSolution[i];
            ConvertedTempSolution[i] = temp;
        }
        
        return ReturnInitialSolutionAsArray(ConvertedTempSolution);
    }

     private static int[] SwapTwoExamsPositionsRandomly (int[] TempSolution, int MaxIterations) 
     {
        
        int[] tmp  = Arrays.copyOf(TempSolution, NumEsami);
        Random rnd = ThreadLocalRandom.current();
        int cont = 0;
       	
	while(cont != MaxIterations) {
            cont++;
            int index1 = rnd.nextInt(NumEsami);
            int index2 = rnd.nextInt(NumEsami);
            if(TempSolution[index1] != TempSolution[index2]) {
                int x = tmp[index1];
                tmp[index1] = tmp[index2];
                tmp[index2] = x;
                if(Constraint_Check.CheckEverything(tmp, CommonStudentMatrix)) {//System.out.println("\n\nSwapped Exam "+ index1+" in time slot " + tmp[index2] + " with Exam "+index2 +" in time slot " + tmp[index1] +"\nAt Iteration " + cont +"\n\n" );
                    return tmp;
                }  
                else {
                    tmp = Arrays.copyOf(TempSolution, NumEsami);
                }
            }
        }
        return TempSolution;
    }
    
    private static boolean CheckChildExistance(int[] Child, int[][] Population, int PopulationSize)
    {
        for(int i=0; i<PopulationSize; i++)
        {
            int counter = 0;
            for(int j=0; j<NumEsami; j++)
            {
                if(Population[i][j] == Child[j])
                    counter++;
            }
            //If counter = NumEsami, then they are identical
            if(counter == NumEsami)
            {
                return true;
            }
                
        }
        
        return false;
    }
    
    
    
    private static boolean CheckConflictInTabuList(int[] CurrentRow, int[][] TabuList, int TabuLength, float Penalty, float BestPenalty)
    {
        //This method returns true if CurrentRow was not in conflict in anything present into TabuList
        for(int i=0; i<TabuLength; i++)
        {
            if(TabuList[i][0]==CurrentRow[1] && TabuList[i][1]==CurrentRow[0] && TabuList[i][2]==CurrentRow[2])
            {
                //We check ASPIRATION CRITERIA
                if(Penalty <= BestPenalty)
                {
                    System.out.println("Aspiration Critera was used!");
                    return true;
                }
                else
                    return false;
            }
        }
        return true;
    }
    
    private static int[][] ShiftRowsMatrix(int[][]TabuList,int TabuLength)
    {
        int StartSlot = TabuList[0][0];
        int EndSlot = TabuList[0][1];
        int Exam = TabuList[0][2];
        
        for(int i=1; i<TabuLength; i++)
        {
            int currentStartSlot = TabuList[i][0];
            int currentEndSlot = TabuList[i][1];
            int currentExam = TabuList[i][2];
            TabuList[i][0] = StartSlot;
            TabuList[i][1] = EndSlot;
            TabuList[i][2] = Exam;
            StartSlot = currentStartSlot;
            EndSlot = currentEndSlot;
            Exam = currentExam;
        }
        return TabuList;
    }
    
    private static int[] SimulatedAnnealingParameter(int NumIterazioni, int[] TempArray, int LengthNeighborhood)
    {
        System.out.println("Starting SIMU");
        int Iteration = 0;
        float alpha = 0.999f;
        int L = 10*LengthNeighborhood;
        ///
        float CurrentMIN = 99999999999f;
        
        Solution = TempArray;
        
        //Set initial Temperature
        
        float TPenalty=0f;
        int[] TSol = new int[NumEsami];
        for(int i=0;i<LengthNeighborhood;i++)
        {
            TSol = GenerateRandomNeighbor(Solution, 3);
            TPenalty += CheckPenaltyFinalSolution(TSol);
        }
        
        TPenalty = TPenalty / LengthNeighborhood;
        
        float T = (Math.abs(TPenalty - CheckPenaltyFinalSolution(Solution)))/0.69f;
        
        int counter = 0;
        
        while(counter < NumIterazioni)
        {
            counter++;
            ///TEST////
            if(Penalty < CurrentMIN)
            {
                CurrentMIN = Penalty;
                //System.out.println("Current penalty is " + CurrentMIN);
            }
            //////////////////////
            
            TempSolution = GenerateRandomNeighbor(Solution, 3);
            TempPenalty = CheckPenaltyFinalSolution(TempSolution);
            if(TempPenalty <= Penalty)
            {
                Solution = TempSolution;
                Penalty = TempPenalty;
            }
            else
            {
               
                float p = (float) Math.exp(-(TempPenalty - Penalty)/T);
                Random rand = new Random();
                float CasualP = rand.nextFloat();
                // Check if we have to swap
                if(CasualP < p)
                {
                    Solution = TempSolution;
                    Penalty = TempPenalty;
                }
            }
            
            //Prepare for next cycle
            Iteration ++;
            if(Iteration%L==0)
            {
               T=T*alpha; 
            }
            
            
        }
        return Solution;
    }
    

    private static float CheckPenaltyFinalSolution(int[] solution)
    {
        float penalty= 0.0f;
        for(int i=0; i<NumEsami; i++)
        {
            for(int j=i+1; j<NumEsami; j++)
            {
                int diff = Math.abs(solution[j]-solution[i]);
                if(diff!=0 && diff<=5)
                {
                    Double d = Math.pow(2.0, 5.0-(double)(diff));
                    penalty += (d*(double)CommonStudentMatrix[i][j])/(double)NumStudenti;
                }
            }
        }
        return penalty;
    }
    
    private static ArrayList<Integer>[] ReconvertSolution(int[] solution)
    {
        ArrayList<Integer>[] TempVar = new ArrayList[NumTimeSlot];
        for (int j = 0; j < TempVar.length; j++)
                TempVar[j] = new ArrayList<>();
        
        for(int i=0; i<NumEsami; i++)
        {
            TempVar[solution[i]].add(i);
        }
        
        return TempVar;
    }
    
    private static void CalculateCommonStudentMatrix()
    {
        CommonStudentMatrix = new int[NumEsami][NumEsami];
        
        for(int i=0; i<NumStudenti; i++)
        {
            for(int j=0; j<NumEsami; j++)
            {
                if(StudentExamMatrix[i][j]==true)
                {
                    for(int k=j; k<NumEsami; k++)
                    {
                        if(StudentExamMatrix[i][k]==true)
                        {
                            CommonStudentMatrix[k][j]++;
                            if(j!=k)
                            {
                                CommonStudentMatrix[j][k]++;
                            }
                        }
                    }
                }
            }
        }
                
    }
    
    private static ArrayList<Exam> SortExamQueue(ArrayList<Exam> ExamQueue, ArrayList<Integer>[] TempSolution)
    {
        //First we count the free slots for each exam
        for(Exam e : ExamQueue)
        {
            int counter = 0;
            e.clearFreeSlots();
            
            for(int index = 0; index < TempSolution.length; index++)
            {
                //IF CELL IS EMPTY
                if(TempSolution[index].isEmpty())
                {
                    counter++;
                    e.addFreeSlot(index);
                }
                //IF CELL IS ALREADY FULL BUT NO CONFLICTS 
                else if(CheckConflictInTempSolution(TempSolution[index],e.getExamID())==true)
                {
                    counter++;
                    e.addFreeSlot(index);
                }
                //IF WE CANNOT PUT THE EXAM INTO THIS CELL
                else
                {
                    //
                }
            }
            //We save the number of possible free slots into the exam
            e.setFreeSlots(counter);  
        }
        
        //Then we sort the exam basing on free slots
        // Sorting
        Collections.sort(ExamQueue, new Comparator<Exam>() {
            @Override
            public int compare(Exam e1, Exam e2)
            {
                return  e1.getFreeSlots().compareTo(e2.getFreeSlots());
            }
        });
        
        //We return
        return ExamQueue;
    }
    
    private static int[] GenerateInitialSolution(int maxIterations, boolean SolutionStore)
    {
        boolean flag = true;    //flag is used to know if a solution was found or not
        ArrayList<Integer>[] TempSolution = new ArrayList[NumTimeSlot];
        ArrayList<Exam> ExamQueue; 
        int iteration;
       
        for(iteration = 0; flag == true && iteration < maxIterations; iteration++)      //Every iteration here starts from scratch
        {
            ExamQueue = new ArrayList<>(Arrays.asList(ExamArray));
            Random randomGenerator = new Random();
            for (int j = 0; j < TempSolution.length; j++)
                TempSolution[j] = new ArrayList<>();
            flag = false;
        
            while(!ExamQueue.isEmpty() && flag == false)
            {
                //We remove the first exam, which will be the one with the least amount of free possible slots
                ExamQueue = SortExamQueue(ExamQueue,TempSolution);
                Exam CurrentExam = ExamQueue.get(0);
                ExamQueue.remove(0);

                //Generation of a random integer between 0 and the number of free slots for the current exam
                int numSlot = CurrentExam.getFreeSlots();

                //if available slots exist, insert the exam in the randomly generated one
                if(numSlot != 0)  {
                    int freeSlotIndex = randomGenerator.nextInt(numSlot);
                    TempSolution[CurrentExam.getSlot(freeSlotIndex)].add(CurrentExam.getExamID());
                }
                else {
                    flag = true;        //Flag becomes true if current iteration is blocked because has nowhere to go (no available free slots for him)
                }
            }
        }
        
        if(flag)
        {
            //No solution was found in maxIterations
            System.out.println("None of the " + maxIterations + " iteration resulted in a valid solution");
            //If method returns -1, that means no initial solution was found in maxIterations
            return null;
        }
        else
        {
            //A solution was found
            //First we calculate its penalty
            CheckPenalty(TempSolution);
            //WE EVENTUALLY CONVERT THE SOLUTION INTO AN EXAM ARRAY storing it into Solution
            if(SolutionStore) {
               ConvertInitialSolution(TempSolution);
            }
            //We return the number of the iteration at which the solution was found
            return ReturnInitialSolutionAsArray(TempSolution);
        }
            
    }
    
    private static boolean CheckConflictInTempSolution(ArrayList<Integer> Cell,int CurrentExam)
    {
        for(int i : Cell)
        {
            if(CommonStudentMatrix[i][CurrentExam]!=0)
                return false;
        }
        return true;
    }
    
    private static void ConvertInitialSolution(ArrayList<Integer>[] TempSolution)
    {        
        for(int i=0; i<TempSolution.length; i++)
        {
            //IF a timeslot is empty
            if(TempSolution[i].isEmpty())
            {
                continue;
            }
            //IF a timeslot is NOT empty
            else
            {
                for(int exam : TempSolution[i])
                {
                    Solution[exam] = i;
                }
            }
        }
    }
    
    private static int[] ReturnInitialSolutionAsArray(ArrayList<Integer>[] TempSolution) 
    {
        
        int[] Solution = new int[NumEsami];
        
        for(int i=0; i<TempSolution.length; i++)
        {
            //IF a timeslot is empty
            if(TempSolution[i].isEmpty())
            {
                continue;
            }
            //IF a timeslot is NOT empty
            else
            {
                for(int exam : TempSolution[i])
                {
                    Solution[exam] = i;
                }
            }
        }
        return Solution;
    }
    
    private static void Calculate_Exam_Array()
    {
        ExamArray = new Exam[NumEsami];
        
        for(int i=0; i<NumEsami; i++)
        {
            int count_conflict = 0;
            int count_penalty = 0;
            for(int j=0; j<NumEsami; j++)
            {
                if(i!=j && CommonStudentMatrix[i][j]!=0)
                {
                    count_conflict++;
                    count_penalty = count_penalty + CommonStudentMatrix[i][j];
                }
            }
            
            ExamArray[i] = new Exam(i,count_conflict,count_penalty);
            
        }
        
    }
    
    private static void CheckPenalty(ArrayList<Integer>[] TempSolution)
    {
        
        float penalty = 0;
        for(int i=0; i<NumTimeSlot; i++)
        {
            for(int exam:TempSolution[i])
            {
                for(int k=i+1; k<i+6 && k<NumTimeSlot; k++)
                {
                    for(int varExam:TempSolution[k])
                    {
                        Double d = Math.pow(2.0, 5.0-(double)(k-i));
                        penalty += (d*(double)CommonStudentMatrix[exam][varExam])/(double)NumStudenti;
                    }
                }
            }
        }
        Penalty = penalty;
    }
    
    private static void ReadFromFile() throws FileNotFoundException, IOException
    {
        
        //LETTURA NUMERO TIME SLOT
        BufferedReader br = new BufferedReader(new FileReader(InstanceName+".slo"));
        try {
            String line = br.readLine();
            NumTimeSlot = Integer.parseInt(line);
        } finally {
            br.close();
        }
        
        //LETTURA NUMERO ESAMI
        br = new BufferedReader(new FileReader(InstanceName+".exm"));
        try {
            String line = br.readLine();
            while (line != null) {
                line = br.readLine();
                NumEsami++;
            }
        } finally {
            br.close();
        }
        
        //LETTURA STUDENTI-ESAMI
        br = new BufferedReader(new FileReader(InstanceName+".stu"));
        try {
            
            ArrayList<boolean[]> TempMatrix = new ArrayList<>(); 
            String OldStudent = null;
            boolean[] TempBoolArray = null;
            int IDesame = -1;
            String line = br.readLine();
            
            while (line != null) {
                
                String[] SplitRes = line.split(" "); //Primo elemetno è studente, secondo id esame
                IDesame = Integer.parseInt(SplitRes[1]);
                //CASO CAMBIO STUDENTE
                if(OldStudent == null)
                {
                    TempBoolArray = new boolean[NumEsami];
                    TempBoolArray[IDesame - 1] = true;
                    NumStudenti++;
                    
                }
                else if(!OldStudent.equals(SplitRes[0]))
                {
                    TempMatrix.add(TempBoolArray);
                    TempBoolArray = new boolean[NumEsami];
                    TempBoolArray[IDesame - 1] = true;
                    NumStudenti++;
                }
                else //STUDENTE VECCHIO = STUDENTE NUOVO
                {
                    TempBoolArray[IDesame - 1] = true;
                }
                
                //Prepare for next iteration
                OldStudent = SplitRes[0];
                line = br.readLine();
            }
            
            //FINAL ITERATION
            TempBoolArray[IDesame - 1] = true;
            TempMatrix.add(TempBoolArray);
            
            //CONVERSIONE DA ARRAYLIST A MATRICE
            StudentExamMatrix = TempMatrix.toArray(new boolean[][] {});
                        
        } finally {
            br.close();
        }
    }
    
}
