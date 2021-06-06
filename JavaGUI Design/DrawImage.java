package project;

import javax.swing.*;
import java.awt.*;
public class DrawImage extends JPanel {
    private Image image;
    public DrawImage(Image im) {
		image=im;
    }
    public void paintComponent(Graphics g) {
    	super.paintComponent(g);
		g.drawImage(image,0,0,getWidth(),getHeight(),Color.white,this);
    }
    public void setImage(Image im) {
		image=im;
		repaint();
    }
    public Image getImage(){
    	return image;
    }
}
