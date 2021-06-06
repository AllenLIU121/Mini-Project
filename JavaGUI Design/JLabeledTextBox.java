package project;

import javax.swing.text.JTextComponent;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

public class JLabeledTextBox extends JPanel{
	
		Vector listText;
		Vector listLabel;
		
		public static final int TEXTFIELD=0;
		public static final int TEXTAREA=1;
		public static final int PASSWORD=2;
	/// Constructeur 1
	public JLabeledTextBox(){
		listText= new Vector();
		listLabel= new Vector();
		//setLayout(new BoxLayout(this,BoxLayout.Y_AXIS));	
		setLayout(new GridLayout(0,1));
	//	setLayout(new BorderLayout());	
	}
	
	///Constructeur 2
	public JLabeledTextBox(String[] Labels){
		
			this();
			for(int i=0;i<Labels.length;i++){
				add(Labels[i]);
			}
	}
	
	///////////////////////////////////////////////////
	public void add(String label){
		
		JLabel l= new JLabel(label);
		listLabel.add(l);
		JTextComponent t= new JTextField();
		listText.add(t);
		
		JPanel pan = new JPanel();
		pan.setLayout(new FlowLayout(FlowLayout.LEFT));
		
		pan.add(l);
		pan.add(t);
		
		add(pan);
	}
	
	public void addTextComponent(String label,int choice){
		JTextComponent t;
		switch(choice){
			case TEXTFIELD : t= new JTextField();break;
			case TEXTAREA  : t= new JTextArea();break;
			case PASSWORD  : t= new JPasswordField();break;
			default : t= new JTextField();break;
		}
		addToPanel(label,t);
	}
	private void addToPanel(String label,JTextComponent t){
		JLabel l= new JLabel(label);
		listLabel.add(l);
		listText.add(t);
		
		JPanel pan = new JPanel();
		pan.setLayout(new FlowLayout(FlowLayout.LEFT));
		
		pan.add(l);
		pan.add(t);
		
		add(pan);
	}
	public void add(String label,int len){
		
		JLabel l= new JLabel(label);
		listLabel.add(l);
		JTextComponent t= new JTextField(len);
		listText.add(t);
		
		JPanel pan = new JPanel();
		pan.setLayout(new FlowLayout(FlowLayout.LEFT));
		
		pan.add(l);
		pan.add(t);
		
		add(pan);
		
		}
		
		
		
		public void add(String label,int len,boolean state){
			
			JTextComponent t= new JTextField(len);
			JLabel l = new JLabel(label);
			listText.add(t);
			listLabel.add(l);
			
			JPanel pan = new JPanel();
			pan.setLayout(new FlowLayout(FlowLayout.LEFT));
			
			pan.add(l);
			pan.add(t);
			
			if(state){
				t.addKeyListener(num);
			}
			else 
				t.removeKeyListener(num);
			
			add(pan); 
		}
		//---------------------------------------
		public void setNum(int Index,boolean state){
			
			JTextComponent t= ((JTextComponent)(listText.get(Index)));
			if(state){
				t.addKeyListener(num);
			}
			else
			t.removeKeyListener(num);

		}
		
		public void setNumAll(boolean state){
			for(int i=0;i<listText.size();i++){
				setNum(i,state);
			}
		}
	///////////////////////////////////////////////////////////////////////////	
		
		public String getValue(int Index){
				
				JTextComponent t = (JTextComponent)(listText.get(Index));
				return t.getText();
				//return ((JTextField)(listText.get(Index))).getText();
		}
	///////////////////////////////////////////////////////////////////////////
	public String getValue(String label){
		int Index=-1;
		for(int i=0;i<listLabel.size();i++){
			JLabel l= (JLabel)(listLabel.get(i));
				if(l.getText().equals(label)){
					
					Index=i;
					break;
				}
		}
		return getValue(Index);
		
	}
	////////////////////////////////////////////////////
	public String[] getValues(){
		String[] s=new String[listText.size()];
		for(int i=0;i<listText.size();i++){
			
		 /*String m=getValue(i);
		 s[i]=m;
		  break;*/
		  s[i]= getValue(i);
		  }
     	return s;
	}
	
	///////////////////////////////////////////////////
	 
	public void setValue(int Index,String Value){
		((JTextComponent)(listText.get(Index))).setText(Value);	
	}
	
	
	public void setValues(String[] Values){
		
		for(int i=0;i<Values.length;i++)
			setValue(i,Values[i]);	
	}
///////////////////////////////////////////////////////
			public static final KeyListener num = new KeyAdapter()	{
			public void keyTyped(KeyEvent e){
			
				if(e.getKeyChar()<'0' || e.getKeyChar()>'9')
					e.setKeyChar((char)0);
			}
		};	

//----------------------------------------------------------------	
		public static final KeyListener Locked = new KeyAdapter(){
		public void keyTyped(KeyEvent e){
			
			e.setKeyChar((char)0);
		}
	};
	
	public void setLocked(int Index,boolean state){
		
		JTextComponent t= ((JTextComponent)(listText.get(Index)));
		if(state){
			t.addKeyListener(Locked);
		}
		else t.removeKeyListener(Locked);
	}
	
	public void setLockedAll(boolean state){
		
			for(int i=0;i<listText.size();i++)
				setLocked(i,state);
	}
	
	public void setEnabled(int Index,boolean state){
		((JTextComponent)(listText.get(Index))).setEnabled(state);
	}
	
	public void setEnabledAll(boolean state){
		for(int i=0;i<listText.size();i++){
			setEnabled(i,state);
		}
	}
	
	/*public void setFocus(int Index){
		JTextField t = new JTextField();
		
	}*/
	public void setLabelPreferredSizeAt(int Index,int height,int width){
		((JLabel)(listLabel.get(Index))).setPreferredSize(new Dimension (height,width));
	}
	
	public void setLabelsPreferredSize(int height,int width){
		for(int i=0;i<listLabel.size();i++){
			setLabelPreferredSizeAt(i,height,width);
		}
	}
	
		public void setTextPreferredSizeAt(int Index,int height,int width){
		((JTextComponent)(listText.get(Index))).setPreferredSize(new Dimension (height,width));
	}
	
	public void setTextsPreferredSize(int height,int width){
		for(int i=0;i<listText.size();i++){
			setTextPreferredSizeAt(i,height,width);
		}
	}
	public void clear(){
		for(int i=0;i<listText.size();i++){
			((JTextComponent)(listText.get(i))).setText("");
		}
	}
	public int FieldsCount(){
		return listText.size();
	}
////////////////////////////////   MAIN  ////////////////////////////////////////////:	
	/*public static void main(String[] args){

		JFrame f= new JFrame();
		f.setBounds(200,100,500,300);
		
	 //	String	s[]={"Nom","Prenom","age"};
		LabelTextBox lbtext = new LabelTextBox();
		
		 lbtext.add("Nom");
		 lbtext.add("Prenom");
		 lbtext.add("Age");
		 lbtext.add("Adresse");
		 
		lbtext.add("Telephone",20,true);
		lbtext.setLabelsPreferredSize(50,50);
		lbtext.setTextsPreferredSize(40,40);
		System.out.println(lbtext.getValue(2));
		lbtext.setLocked(3,true);
		lbtext.setEnabled(1,false);
	
	
	//	JPanel p = new JPanel();
	//	p.setLayout(new BoxLayout(p,BoxLayout.Y_AXIS));
		
	//	p.add("North",lbtext);
	
		f.add(lbtext);
		f.setVisible(true);
	}*/
	
}