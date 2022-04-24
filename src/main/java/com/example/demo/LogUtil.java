package com.example.demo;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

public class LogUtil {
    private BufferedWriter writer;

    public LogUtil(){
        try{
            String fileName = "dsv_log_"+LocalDateTime.now().format(DateTimeFormatter.ofPattern("dd-M-yyyy HH-mm-ss")) + ".txt";
            this.writer = new BufferedWriter(new FileWriter(fileName, false));
        }catch(IOException ioe){
            System.out.println(ioe.getMessage());
        }

    }

    public void log(String message){
        try{
            writer.write(message + "\n");
        }catch (IOException ioe){
            System.out.println(ioe.getMessage());
        }
    }

    public void close(){
        try{
            writer.close();
        }catch (IOException ioe){
            System.out.println(ioe.getMessage());
        }
    }

}
