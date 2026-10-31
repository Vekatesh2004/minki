# DNA Helix Background Animation Feature

## Overview
Added a pure CSS animated DNA double helix background to enhance the visual appeal of the HelixAI platform without impacting performance.

## Implementation Date
August 2, 2026

## Technical Details

### What Was Added
1. **CSS Classes** (lines ~828-850 in `simple_backend.py`)
   - `.dna-helix` - Container for the animation
   - `.dna-strand` - Individual DNA strands (3 strands at different positions)
   - `.dna-base-pair` - Individual base pairs with rotating animation
   - `.dna-base` - Base pair spheres with gradient colors
   - `.dna-connector` - Connecting lines between base pairs

2. **Animations**
   - `rotate3d` - 8-second 3D rotation for each base pair
   - `dnaFloat` - 20-second floating effect for entire strands
   - Staggered delays for natural helical appearance

3. **HTML Structure** (after line 1041)
   - Three DNA strands positioned at 15%, 50%, and 85% horizontally
   - 12 base pairs per strand
   - Each base pair has two bases connected by a line

### Visual Characteristics
- **Opacity**: 8% (very subtle, doesn't distract)
- **Colors**: Gradient from blue (#2d68df) → teal (#14a19a) → purple (#765bd9)
- **Position**: Fixed background, behind all content (z-index: -1)
- **Performance**: Pure CSS, no JavaScript, minimal CPU usage

### Accessibility
- Properly labeled with `aria-hidden="true"`
- Respects `prefers-reduced-motion` - animation is disabled for users who prefer reduced motion
- Does not interfere with screen readers

### Browser Compatibility
- Works on all modern browsers (Chrome, Firefox, Safari, Edge)
- Gracefully degrades on older browsers
- Uses standard CSS transforms and animations

## Why Pure CSS Instead of Video?
- **Performance**: No file downloads, instant load
- **Size**: ~2KB CSS vs 100+ MB for 4K video
- **Compatibility**: Works everywhere, no codec issues
- **Maintainability**: Easy to customize colors/speed
- **Copyright**: No licensing concerns

## Testing
Test on local server: http://127.0.0.1:8000
1. Open the page - you should see subtle rotating DNA helices in the background
2. Scroll and interact - animation continues smoothly
3. Test with example file or upload - functionality unchanged
4. Check browser console - no errors

## Customization Options
To adjust the animation, modify in `simple_backend.py`:

### Make it more/less visible
```css
.dna-helix { opacity: .08; }  /* Change from .08 to .15 for more visible */
```

### Change speed
```css
@keyframes rotate3d { ... }  /* Change from 8s to 12s for slower */
@keyframes dnaFloat { ... }  /* Change from 20s to 30s for slower */
```

### Change colors
```css
.dna-base { 
    background: linear-gradient(135deg, #2d68df, #14a19a); /* Modify colors */
}
```

### Disable completely
Remove the `.dna-helix` div from HTML (around line 1070)

## Files Modified
- `simple_backend.py` - Added CSS classes and HTML structure

## Deployment
After testing locally:
1. Transfer file to production server
2. Restart server: `sudo systemctl restart pharmacogenomics`
3. Hard refresh browser: `Ctrl + Shift + R`

## Notes
- Original SVG background (`science-bg`) is still present and works in harmony
- Both backgrounds create a layered depth effect
- Disease annotation and disease badges features are preserved
- No impact on existing functionality
