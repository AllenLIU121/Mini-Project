package project;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.util.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.awt.*;
public class JToolBarButton extends JToolBar{
		
	private	Vector TB;
	
	public JToolBarButton(String[] labels){
		TB= new Vector();
		for(int i=0;i<labels.length;i++){
				ImageIcon ig = new ImageIcon(new ImageIcon("Images/"+
				labels[i]+".png").getImage().getScaledInstance(50, 50, Image.SCALE_AREA_AVERAGING));
				
				JButton b= new JButton(ig);
				b.setToolTipText(labels[i]);
				b.setName(labels[i]);
				//b.setText(labels[i]);
				TB.add(b);
				add(b);
		}
	}
	
	public JToolBarButton(String [] icons,String[] ToolTips){
	
		TB= new Vector();
		for(int i=0; i<icons.length;i++){
			
			JButton b = new JButton((new ImageIcon("images/"+icons[i])));
			b.setToolTipText(ToolTips[i]);
			b.setName(ToolTips[i]);
			TB.add(b);
			add(b);
		}
	
	}
	
	public void enabledLabels(boolean state){
			for(int i=0; i<TB.size();i++){
					JButton b= (JButton)TB.get(i);
				if(state){
				
					b.setText(b.getToolTipText());	
				}
				else 
					b.setText(null);
				}
		
	}
	
	public void setToolTips(boolean state){
		for(int i=0;i<TB.size();i++){
				JButton b = (JButton)TB.get(i);
			if(state){
			
				b.setToolTipText(b.getText());
			}
			
			else
				b.setToolTipText(null);
		}
	}
	
	
	
	
	public void addActionListener(ActionListener al){
		for (int i=0;i<TB.size();i++){
			((JButton)TB.get(i)).addActionListener(al);
		}
	}
	
	
}