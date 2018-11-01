/* $Id$ $Revision$ */
/* vim:set shiftwidth=4 ts=8: */

/*************************************************************************
 * Copyright (c) 2011 AT&T Intellectual Property 
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors: See CVS logs. Details at http://www.graphviz.org/
 *************************************************************************/

#import <Foundation/Foundation.h>
#import <AppKit/AppKit.h>
#import <Quartz/Quartz.h>

@interface GVWindowController : NSWindowController <NSUserInterfaceValidations>
{
	IBOutlet PDFView *documentView;
}

- (id)init;
- (void)setDocument: (NSDocument *)document;
- (void)awakeFromNib;

- (void)graphDocumentDidChange:(NSNotification*)notification;
- (NSRect)windowWillUseStandardFrame:(NSWindow *)window defaultFrame:(NSRect)defaultFrame;

- (IBAction)actualSizeView:(id)sender;
- (IBAction)zoomInView:(id)sender;
- (IBAction)zoomOutView:(id)sender;
- (IBAction)zoomToFitView:(id)sender;

- (IBAction)printGraphDocument:(id)sender;

- (void)PDFViewWillClickOnLink:(PDFView *)sender withURL:(NSURL *)URL;

- (BOOL)validateUserInterfaceItem:(id <NSValidatedUserInterfaceItem>)anItem;

- (void)dealloc;
@end
