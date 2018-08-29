/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package oma_assignment;

import java.util.ArrayList;

/**
 *
 * @author Arturo
 */
public class Exam {

    private int ExamID,NumConflicts,NumCommonStudents,FreeSlots;
    private ArrayList<Integer> AvailableSlots;
    
    public Exam(int examID, int numConflicts, int numCommonStudents)
    {
        ExamID = examID;
        NumConflicts = numConflicts;
        NumCommonStudents = numCommonStudents;
        FreeSlots = -1;
        AvailableSlots = new ArrayList();
    }
   
    ///////// ARRAYLIST UTILITIES /////////
    
    public void clearFreeSlots() {
        AvailableSlots.clear();
    }
            
    public void addFreeSlot(int slot) {
        AvailableSlots.add(slot);
    }
    
    public ArrayList<Integer> returnSlots()
    {
        return AvailableSlots;
    }
    
    public void printFreeSlots() {
        AvailableSlots.stream().forEach(Element -> System.out.print(Element + " "));
        System.out.println();
    }
    
    public int getSlot(int i) {
        return AvailableSlots.get(i);
    }
    
    //////////// NORMAL GETTERS & SETTERS ///////////////

    public Integer getFreeSlots() {
        return AvailableSlots.size();
    }
    
    public void setFreeSlots(int FreeSlots) {
        this.FreeSlots = FreeSlots;
    }

    public int getExamID() {
        return ExamID;
    }

    public int getNumConflicts() {
        return NumConflicts;
    }

    public int getNumCommonStudents() {
        return NumCommonStudents;
    }
  
}
