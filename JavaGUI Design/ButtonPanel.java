package project;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*; 



public class  ButtonPanel extends JPanel{
	private Vector ButtonList;
	
	public static final int FlowCenter = 0;
	public static final int FlowRight  = 1;
	public static final int FlowLeft   = 2;
	public static final int BoxYAXIS   = 3;
	public static final int BoxXAXIS   = 4;
	
	
	public ButtonPanel(String labels[],int type){
		super();
		
		ButtonList=new Vector();
	
		switch(type){
			case FlowCenter : setLayout(new FlowLayout(FlowLayout.CENTER));break;
			case FlowRight : setLayout(new FlowLayout(FlowLayout.RIGHT));break;
			case FlowLeft : setLayout(new FlowLayout(FlowLayout.LEFT));break;
			case BoxYAXIS : setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));break;
			case BoxXAXIS : setLayout(new BoxLayout(this, BoxLayout.X_AXIS));break;		
			}
			
		createButton(labels);
	}
	public ButtonPanel(String[] labels,int row,int col){
		super();
		ButtonList=new Vector();
		setLayout(new GridLayout(row,col));
		createButton(labels);
		
	}
	
	public ButtonPanel(String[] labels,int row,int col,int Xspace,int Yspace){
		super();
		ButtonList=new Vector();
		setLayout(new GridLayout(row,col,Xspace,Yspace));
		createButton(labels);
		
	}
	
	private void createButton(String[] labels){
		for(int i=0;i<labels.length;i++){
			JButton b= new JButton(labels[i]);
			b.setName(labels[i]);
			ButtonList.add(b);
			add(b);		
		}
		
	}
	
	public void setPreferdSizeAt(int Index,int Height,int Width){
	
		((JButton)ButtonList.get(Index)).setPreferredSize(new Dimension(Height, Width));
	}
	
	
	public void setPreferdSize(int Height,int Width){
		
		for(int i=0;i<ButtonList.size();i++)
			setPreferdSizeAt(i, Height,Width);
			
	}
	
	
//////////////////////////////////////////////////
	
	public void addActionListener(ActionListener al){
		for(int i=0;i<ButtonList.size();i++){
			((JButton)ButtonList.get(i)).addActionListener(al);
			}
	}
	
	
	
}