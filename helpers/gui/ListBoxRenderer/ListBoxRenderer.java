import java.awt.*;
import javax.swing.*;
public class ListBoxRenderer extends JLabel implements ListCellRenderer
{
    // Class constructors
    public ListBoxRenderer() {
        setOpaque(true);
        setHorizontalAlignment(LEFT);
        setVerticalAlignment(CENTER);
    }

    // Return label and set tooltip
    public Component getListCellRendererComponent(
            JList list,
            Object value,
            int index,
            boolean isSelected,
            boolean cellHasFocus)
    {
        if (isSelected) {
            // Selected cell item
            setBackground(list.getSelectionBackground());
            setForeground(list.getSelectionForeground());
        } else {
            // Unselected cell item
            setBackground(list.getBackground());
            setForeground(list.getForeground());
        }

        String label = value.toString();
        list.setToolTipText(label);
        setText(label);

        return this;
    }
}
