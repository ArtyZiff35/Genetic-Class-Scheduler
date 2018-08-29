/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package oma_assignment;

/**
 * @author Arturo
 */
public class Constraint_Check {
    
    public static boolean CheckEverything(int[] Solution,int[][] CommonStudentMatrix)
    {
        int[][] DistanceMatrix = CalculateDistanceMatrix(Solution);
        if(!CheckOverlapping(DistanceMatrix,CommonStudentMatrix))
            return false;
        
        //To check last constraint, we need to see if a time slot is put twice into the same array slot
        
        return true;
    }
    
    public static boolean CheckOverlapping(int[][] DistanceMatrix,int[][] CommonStudentMatrix)
    {
        for(int i=0; i<DistanceMatrix[0].length; i++)
        {
            for(int j=i+1; j<DistanceMatrix[0].length; j++)
            {
                if(DistanceMatrix[i][j]==0)
                    if(CommonStudentMatrix[i][j] != 0)
                        return false;
            }
        }
        return true;
    }
   
    public static int[][] CalculateDistanceMatrix(int[] Solution)
    {
        int Dim = Solution.length;
        int[][] TempMatrix = new int[Dim][Dim];
        
        for(int i=0; i<Dim; i++)
        {
            for(int j=i+1; j<Dim; j++)
            {
                int Diff = Solution[i]-Solution[j];
                TempMatrix[i][j] = Math.abs(Diff);
                TempMatrix[j][i] = Math.abs(Diff);
            }
        }
        
        return TempMatrix;
    }
    
}
