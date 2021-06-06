package project;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JPanel;

public class JButtonImage extends JPanel{
	
	private Vector ButtonList;
	
	public static final int FlowCenter = 0;
	public static final int FlowRight  = 1;
	public static final int FlowLeft   = 2;
	public static final int BoxYAXIS   = 3;
	public static final int BoxXAXIS   = 4;
	
	private JButtonImage(String []icones,String []Tooltips,int choice){
		ButtonList=new Vector();
	
		switch(choice){
			case FlowCenter : setLayout(new FlowLayout(FlowLayout.CENTER));
			case FlowRight : setLayout(new FlowLayout(FlowLayout.RIGHT));
			case FlowLeft : setLayout(new FlowLayout(FlowLayout.LEFT));
			case BoxYAXIS : setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
			case BoxXAXIS : setLayout(new BoxLayout(this, BoxLayout.X_AXIS));		
		}
			
		createButton(icones,Tooltips);
	}
	
	public JButtonImage(String []icones,int choice){
		this(icones,icones,choice);
	}
	private void createButton(String[] icones,String []Tooltips){
		for(int i=0;i<icones.length;i++){
			JButton b= new JButton(new ImageIcon("Images/"+icones[i]+".gif"));
			b.setName(Tooltips[i]);
			b.setToolTipText(Tooltips[i]);
			add(b);
			ButtonList.add(b);
		}
	}
	
	public void setPreferdSizeAt(int Index,int Height,int Width){
		((JButton)ButtonList.get(Index)).setPreferredSize(new Dimension(Height, Width));
	}
	
	
	public void setPreferdSize(int Height,int Width){
		
		for(int i=0;i<ButtonList.size();i++)
			setPreferdSizeAt(i, Height,Width);
			
	}
	
	public JButtonImage(String[] icones,String []Tooltips,int row,int col){
		super();
		ButtonList=new Vector();
		setLayout(new GridLayout(row,col));
		createButton(icones,Tooltips);
		
	}
	
	public JButtonImage(String[] icones,int row,int col){
		this(icones,icones,row,col);
	}
	
	public JButtonImage(String[] icones,String []Tooltips,int row,int col,int Xspace,int Yspace){
		super();
		ButtonList=new Vector();
		setLayout(new GridLayout(row,col,Xspace,Yspace));
		createButton(icones,Tooltips);
	}
	public void addActionListener(ActionListener al){
		for(int i=0;i<ButtonList.size();i++){
			((JButton)ButtonList.get(i)).addActionListener(al);
		}
	}
		
}