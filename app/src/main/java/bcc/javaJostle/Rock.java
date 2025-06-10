package bcc.javaJostle;

import java.util.ArrayList;

public class Rock extends Robot {
    public Rock(int x, int y){
        super(x, y, 5, 1, 3, 1,"Rock", "rock.png", "rock.png");
        
        // Health: 3, Speed: 2, Attack Speed: 2, Projectile Strength: 3
        // Total = 10
    }

    public void think(ArrayList<Robot> robots, ArrayList<Projectile> projectiles, Map map, ArrayList<PowerUp> powerups) {
        //rock robot is not smart and doesn't think very well. 
                
    }

}
